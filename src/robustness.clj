(ns robustness
 (:require [clojure.contrib.string :as str]
           [clojure.java.shell :as shell]
           [clojure.contrib.io :as io]
           [incanter.stats :as stats]
           [incanter.charts :as charts]
           [clojure.contrib.json :as json]
           )
  (:use edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.snippets-math
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.core :only (view)]
        refold
        ))

(defn fold
  "Folds a sequence of RNA and returns only the target structure."
  
  [s]
  (->> ((shell/sh "RNAfold"
                  "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                  "--noPS"
                  :in s)
        :out)
       (str/split-lines)
       second
       (str/split #" ")
       first))

(defn mutant-neighbor
  "Takes a string s and finds the 3L 1-mer mutants. String can only contain
   letters A, C, G, U."
  
  [s]
  (let [s (.toUpperCase s)]
    (flatten
     (for [i (range (count s))]
       (map (fn [r]
              ;;makes new sequence with the substitution r
              (str (subs s 0 i) r (subs s (inc i) (count s))))
            (keys (dissoc {"A" 1 "G" 1 "U" 1 "C" 1} ;;3 other bases to sub
                          (subs s i (inc i)))))))))

(defn inverse-fold
  "Given a target structure, it will use RNAinverse to find n sequences which fold into an identical structure"
  
  [target n]
  (loop [c 0
         cand []]
    (if (< c n)
      (let [x (remove nil?
                      (flatten 
                       (map (fn [[s ensemble]]
                              (when-not (re-find #"d=" s) (re-find #"\w+" s)))
                            (->> ((shell/sh "RNAinverse"
                                            "-Fmp"
                                            (str "-R" n)
                                            "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                            :in target)
                                  :out)
                                 (str/split-lines)
                                 (partition-all 2))
                            )))]
        (recur (count (distinct cand))
               (concat (distinct cand) x)))
      (take n (distinct cand)))))

(defn suboptimals
  "Finds the centroid structure of suboptimal structures and a vector
   representation using 0's and 1's using a RNAmutants or RNAsubopt. s
   is the RNA sequence (case doesn't matter as it will be all
   upper-cased) and n is the number of suboptimal structures to
   consider."
  
  [s n]
  (let [;;s "AACGAUCCCGGG"
        ;;n 10
        s (.toUpperCase s)
        RNAmutants 0 #_((shell/sh "./RNAmutants"
                              "-l" "./lib/"
                              "--mutation" "1"
                              "-n" (str n)
                              "--input-string" s
                              :dir "/home/kitia/Desktop/RNAmutants/")
                    :out)
        RNAsubopt ((shell/sh "RNAsubopt"
                             "-p" (str n) ;samples according to
                                          ;Boltzmann distribution
                             :in s)
                   :out)
        out (->> RNAsubopt str/split-lines)
        structures (->> (if (some #(re-find #"\w+" %) out)
                          (drop-until #(re-find #"\> sampling \d+" out))
                          out)
                        (remove #(re-find #"[^\(\)\.]" %))
                        )
        struct->matrix (fn [st] (reduce (fn [m kv]
                                         (assoc m kv 1))
                                       {} (make_pair_table st)))
        Z->centroid (fn [matrix]
                      (sort-by key
                               (reduce (fn [m [[i j] p]]
                                         (assoc m i "(" j ")"))
                                       (into {}
                                             (map #(vector %1 %2) (range (count s)) (repeat ".")))
                                       (filter (fn [[[i j] p]]
                                                 (and (< i j)
                                                      (>= p 0.5)))
                                               matrix))))
        struct->vector (fn [st] (map #(if (= \. %) 0 1) (seq st)))
        partition-function (reduce  (fn [m [k v]]
                                      (assoc m k (/ v n)))
                                    {}  (apply merge-with + (map struct->matrix structures)))
        centroid (apply str (vals (Z->centroid partition-function)))]
    ;;(doseq [i structures] (prn i))
    ;;(prn (apply str (vals (Z->centroid partition-function))))
    [centroid (struct->vector centroid)]
    ))

(defn neutrality
  "takes a string s and  returns neutrality of the seq  when compared
   to each of the 3L 1-mutant neighbors"

  [s & [{foldtype :foldtype cons :cons :or {foldtype "mfe" cons nil}}]]
  (let [st (cond
            (= foldtype "mfe")
            (fold s)

            (= foldtype "centroid")
            (first (suboptimals s 10000)))
        neighbors (mutant-neighbor s)
        neutfn (fn [x] (/ (- (count s) x)
                       (count s)))
        ]
    (if-let [x cons]
      (fn [cons]
        (let [n 7
              cons (->> ((profile (read-sto cons)) :structure)
                        first
                        (str/replace-re #"\:|\-" "." ))]
          (map (fn [neighbor]
                 (let [js (jensen-shannon (probs n cons) (probs n (fold neighbor)))
                       norm (levenshtein st (fold neighbor))
                       c (apply levenshtein (-> (sw cons (fold neighbor))
                                                first
                                                rest))
                       ]
                   [(neutfn norm) (neutfn c) (- 1 js)]))
               neighbors)))
      (fn []
        (map (fn [neighbor]
               (neutfn (levenshtein st (fold neighbor))))
             neighbors)
        )
      )))

(defn neutrality_rand
  "neutrality from finding a random seq with the same structure. if neutrality
   for native seq > rand seq then the native seq is more robust."
  
  ([target]
     (neutrality_rand target 100))
  
  ([target n]
     (let [cand (inverse-fold target n)]
       #_(prn distinct? cand)
       ;;(filter #(= target (fold %)) cand)
       (map #(stats/mean (neutrality %)) cand))))




(defn psdc [dataset]
  (map (fn [wt mut]
         (map #(* (Math/sqrt (count wt))
                  (- 1 (pearsonsCC wt %))) mut))
       (first dataset) (second dataset)))

(defn generate-vectors
  "Generates suboptimal structures to find the centroid
   structure (structure and vector representation) of sequence s. n
   inverse folded sequences will be generated based on the centroid
   structure. Then the centroid structure for all 3L 1-mutant
   neighbors will be found for s and its n inverse folded sequences.
   Returns the structure of s then the mutuant in vector form."

  [s n]
  (let [s (.toUpperCase s)
        [stc stv] (suboptimals s 10000)
        wts  (inverse-fold stc n)
        muts (pxmap
              (fn [i] (doall (map #(second (suboptimals % 10000)) i)))
              10
              (map #(mutant-neighbor %) (cons s wts)))] 
    [(cons stv (map #(second (suboptimals % 10000)) wts)) muts]))


(defn makechart
  "creates a JSON datastructure for a chart from a map. The JSON
   output file in .js format can then be used to create a highcharts
   chart

   f is typically \"gaisr/robustness/highchart-test.js\""
  
  [f dataset title subtitle]
  (let [;;abc (map #(map (fn [x] (stats/mean x)) (partition-all 3 %))
   ;;barr)
        abc (map #(map (fn [x] (stats/mean x)) (partition-all 3 %)) dataset)
        wt (first abc)
        n (count wt)
        xy {:chart {
                    :renderTo "container",
                    :type "line",
                    :marginRight 130,
                    :marginBottom 50
                    :height 500}
            :title {:text title,
                    :x -20
                    }
            :subtitle {:text subtitle 
                       :x -20}
            :xAxis {:title {:text "position"} 
                    :categories (vec (range 1 (inc n)))}
            :yAxis {:title {:text "pSDC"} 
                    :min 0 :maxPadding 0.001
                    :plotLines [{:value 0
                                 :width 1
                                 :color "#0808080"}]}
            :legend {:layout "vertical"
                     :align "right"
                     :verticalAlign "top"
                     :x -10
                     :y 100
                     :borderWidth 0}
            :series (vec
                     (conj
                      (map (fn [i y]
                             {:name (str "neg" i)
                              :data (vec y)
                              :visible true})
                           (iterate inc 1) (rest abc))
                      {:name "wt"
                       :data (vec wt)}
                      ))}
        ]
    (clojure.contrib.io/with-out-writer f 
      (print "var foo=")
      (prn (json/json-str xy)))))


(defn combinedset []
  (let [avgpt (fn [dataset] (map #(map (fn [x] (stats/mean x)) (partition-all 3 %)) dataset))
        makemap (fn [names data] 
                  (into {}
                        (map (fn [nm x] 0
                               [nm {:wt (first x) :con (rest x)}])
                             names data)))
        exp1 (map avgpt 
                  (map psdc 
                       (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con.clj"))))
        exp1 (makemap [:L10seq2 :L13 :L20 :L21] exp1)
        exp2 (map avgpt
                  (map psdc 
                       (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con2.clj"))))
        exp2 (makemap [:L10seq1 :L10seq2 :L13 :L20 :L21] exp2)                             
        exp3 (map avgpt (read-string (slurp "/home/peis/bin/gaisr/robustness/inverse-struct-con3.clj")))
        exp3 (makemap [:L10seq1 :L10seq2 :L13 :L20 :L21 :FMN] exp3)
        comboexp (merge-with (fn [a b]
                               (let [{curwt :wt curcon :con} a
                                     {newwt :wt newcon :con} b]
                                 (assoc {} :wt (map + curwt newwt) :con (concat curcon newcon))))
                             exp1 exp2 exp3)]
    (reduce (fn [m [k {wt :wt con :con}]]
              (let [n {:L10seq1 2 :L10seq2 3 :L13 3 :L20 3 :L21 3 :FMN 1}
                    w (map (fn [x] (/ x (n k))) wt)
                    c (map stats/mean con)]
                (assoc m k 
                       (assoc {} :wt w :wtavg (stats/mean w) :conavg c
                              :con con))))
            {} comboexp)))

(defn evaluate-inverse-cons4
  "let [foo
   \"/home/peis/bin/gaisr/robustness/inverse-struct-con4.clj\"] Then
   we are evaluating foo which is a vector of maps where each map only
   has 1 key (:L10 :L13 :L20 :L21 :FMN). The pSDC has already been
   found for the wt and each of its 1-mutant neighbors. These values
   still have to be averaged. The ranking of the WT pSDC value shows
   the significance of the pSDC value as compared to 100 inverse-fold
   sequences."

  []
  
  (map (fn [k element]
         (map (fn [x] 
                (let [[wt & muts] (map stats/mean x)] ;averaging
                                                      ;individual pSDC values for 1 sequence
                  [(count (remove #(< wt %) muts)) ;remove mutant pSDC values greater than WT. 
                   (count muts)])) 
              (get element k))) ;multiple sequences for each key. so
                                ;they are mapped over
       (for [i foo] (first (keys i))) foo))

(timefn (fn [] 
                      (let [fdir "/home/peis/bin/gaisr/trainset2/"
                            fs (str/split-lines ((shell/sh "ls" :dir (str fdir "pos/")) :out))
                            fsto fs
                            ffasta (map #(str (str/butlast 3 %) "fasta") fs)
                            x (fn [s cons n] (let [cons (->> ((profile (read-sto cons)) :structure)
                                                             first
                                                             (str/replace-re #"\:|\-" "." ))
                                                   foldtype "centroid"
                                                   st (cond
                                                       (= foldtype "mfe")
                                                       (fold s)
                                                       
                                                       (= foldtype "centroid")
                                                       (first (suboptimals s 10000)))
                                                   neighbors (mutant-neighbor s) ] (prn cons)
                                                   (map (fn [neighbor]
                                                          (let [norm (levenshtein st (fold neighbor))
                                                                c (apply levenshtein (-> (sw cons (fold neighbor))
                                                                                         first
                                                                                         (rest)))
                                                                js (jensen-shannon (probs n cons) (probs n (fold neighbor)))
                                                                neut (fn [x] (/ (- (count s) x)
                                                                                (count s)))]
                                                            [(neut norm) (neut c) (- 1 js)]))
                                                        neighbors)))
                            ]
                        (doall 
                         (map (fn [sto fasta]
                                (let [l (->> (io/read-lines (str fdir fasta))
                                             rest
                                             (take-nth 2)
                                             )]
                                  (for [s l]
                                    (pmap (fn [i]
                                            (let [tx (transpose (x s (str fdir "pos/" sto) i))]
                                              [(pearsonsCC (first tx) (second tx))
                                               (pearsonsCC (first tx) (third tx))
                                               (pearsonsCC (second tx) (third tx))]))
                                          (range 1 11)))
                                  ))
                              fsto ffasta)))))

[((([0.8069075749148139 0.4045922999595296 0.38072041406130336] [0.8069075749148139 0.5274598722165978 0.5027823251365909] [0.8069075749148139 0.5259920891232005 0.5430261861153335]) ([-0.31616569683849693 0.12359703836634316 0.34454319403559786] [-0.31616569683849693 0.09129259901407907 0.37310596448809635] [-0.31616569683849693 0.0018627835837926171 0.3708567503749153]) ([0.8454671331927014 0.5552016610235578 0.568433704535805] [0.8454671331927014 0.656695371035758 0.6144911833195077] [0.8454671331927014 0.7161371601925596 0.6416332718621803]) ([0.9177778900937895 0.5159671430466998 0.6714903261287056] [0.9177778900937895 0.6483723478087349 0.7674213145693864] [0.9177778900937895 0.6862765687488515 0.7984445531044194]) ([0.5377032069919034 0.48734625862535275 0.4730989745839409] [0.5377032069919034 0.3949431421480542 0.25518504778812856] [0.5377032069919034 0.41798408764134093 0.2220047593481736]) ([0.6152842532254943 0.26949996795657455 0.5341083866222135] [0.6152842532254943 0.3769380300378776 0.666065157138977] [0.6152842532254943 0.4576702992273768 0.7214465501755607])) (([0.9171827408433494 0.5114413125455792 0.3606871400640698] [0.9171827408433494 0.46870417138917897 0.316590396979833] [0.9171827408433494 0.4912119644674476 0.32321960664514476]) ([0.020170457736599506 0.5702958130472919 0.08466482208733643] [0.020170457736599506 0.5253672603381837 0.07396393349349417] [0.020170457736599506 0.4196035537933955 0.10591357883103165]) ([0.0960301793527424 0.4388100614207224 0.7393881545642706] [0.0960301793527424 0.40209351779393177 0.7654214563681768] [0.0960301793527424 0.2859668983292572 0.7773992529669893]) ([0.8816558158103343 0.7818932977724412 0.7835324591290789] [0.8816558158103343 0.7756612434337495 0.819313008183512] [0.8816558158103343 0.7551155678050233 0.8203429238613659]) ([-0.7711837145681115 -0.6571250774342878 0.8668784632507986] [-0.7711837145681115 -0.7217468052760775 0.917371260919937] [-0.7711837145681115 -0.773637840769526 0.9557894901982612])) (([0.3856611458673418 0.11355248377046019 -0.157155640773158] [0.3856611458673418 0.13313913904935507 -0.09213030474192538] [0.3856611458673418 0.16080485247064588 -0.09331379065654921]) ([0.9550757236596278 0.6000106271484096 0.756469587528386] [0.9550757236596278 0.6161812684170661 0.7723062158656114] [0.9550757236596278 0.6172242050198491 0.7759496173009961]) ([0.8641601719178124 -0.038082166266373896 0.11231444800985868] [0.8641601719178124 -0.01525110257324259 0.12130583998720473] [0.8641601719178124 0.021672329311787253 0.1156927397043733])) (([0.4373364415168213 -0.3068079823448847 0.018470555831987562] [0.4373364415168213 -0.46119420167835307 0.017167293774204048] [0.4373364415168213 -0.4294303921823325 0.09234497682966229]) ([0.8905408354686025 0.26456554768282237 0.5635466363571816] [0.8905408354686025 0.4482802822543055 0.6929315712438356] [0.8905408354686025 0.5056260369728593 0.7405680844886263]) ([-0.3163374525693546 0.1718190192538901 -0.06993224432009398] [-0.3163374525693546 0.43642593579257116 -0.1125855654929385] [-0.3163374525693546 0.2443040041079991 -0.08883024152577472]) ([0.9070192572463573 0.6328609707458505 0.6701490668915623] [0.9070192572463573 0.4767253780778724 0.5924474666971496] [0.9070192572463573 0.48537616886883567 0.6301285372132759]) ([0.7913241775286965 0.34347985727758 0.5069635782483939] [0.7913241775286965 0.4003439955389045 0.5227337069307542] [0.7913241775286965 0.7195205911952929 0.7326128734326607]) ([-0.2706961177929172 -0.4867060999249683 0.46254864301939] [-0.2706961177929172 -0.47754333433236523 0.4263844779984547] [-0.2706961177929172 -0.31027903339208357 0.3319888018112472]))) 1458451.210406]
