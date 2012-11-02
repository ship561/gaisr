(ns robustness
 (:require [clojure.contrib.string :as str]
           [clojure.java.shell :as shell]
           [clojure.contrib.io :as io]
           [incanter.stats :as stats]
           [incanter.charts :as charts]
           [clojure.contrib.json :as json]
           [clojure.set :as sets]
           [edu.bc.fs :as fs]
           )
  (:use edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.snippets-files
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.utils.snippets-math
        edu.bc.bio.sequtils.dists
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.core :only (view)]
        refold
        smith-waterman
        [consensus_seq :only (profile read-sto change-parens)]
        ))

(def homedir (edu.bc.fs/homedir))

(defn mutant-neighbor
  "Takes a string s and finds the 3L 1-mer mutants. String can only
   contain letters A, C, G, U. with-names returns a vector with the
   mutant name as well"
  
  [s & {:keys [with-names]
        :or {with-names false}}]
  (let [s (.toUpperCase s)]
    (apply concat
           (for [i (range (count s))]
             (map (fn [r]
                    ;;makes new sequence with the base substitution r
                    (if with-names
                      ;;returns the mutation name and seq otherwise just the seq
                      [(str (subs s i (inc i)) i r)
                       (str (subs s 0 i) r (subs s (inc i) (count s)))]
                      (str (subs s 0 i) r (subs s (inc i) (count s)))))
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

(defn struct->matrix
  "creates array of bp locations. Array resembles a hash-map where the
  keys are base-pair locations and the value is 1 if
  present. Locations not present are not represented in the hash-map."
  
  [st]
  (reduce (fn [m kv] ;creates array of bp locations
            (assoc m kv 1))
          {} (make_pair_table st)))

(defn suboptimals
  "Finds the centroid structure of suboptimal structures and a vector
   representation using 0's and 1's using a RNAmutants or RNAsubopt. s
   is the RNA sequence (case doesn't matter as it will be all
   upper-cased) and n is the number of suboptimal structures to
   consider."
  
  [s n & {:keys [centroid-only]
          :or {centroid-only true}}]
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
        Z->centroid (fn [matrix] ;converts the list of subopt
                                ;structures into a centroid
                      (sort-by key
                               (reduce (fn [m [[i j] p]]
                                         (assoc m i "(" j ")"))
                                       (into {}
                                             (map #(vector %1 %2) (range (count s)) (repeat ".")))
                                       ;;keeps bases that have over 50%
                                       ;;representation in suboptimal structures
                                       (filter (fn [[[i j] p]] 
                                                 (and (< i j)
                                                      (>= p 0.5)))
                                               matrix))))
        struct->vector (fn [st] (map #(if (= \. %) 0 1) (seq st))) ;change structure representation
                                                                  ;to a vector of 0's and 1's
        map-structures (map struct->matrix structures) ;sparse matrix
                                                       ;of each structure
        partition-function (reduce  (fn [m [k v]]
                                      (assoc m k (/ v n)))
                                    {}  (apply merge-with + map-structures))
        centroid (apply str (vals (Z->centroid partition-function)))]
    ;;(doseq [i structures] (prn i))
    ;;(prn (apply str (vals (Z->centroid partition-function))))
  (if centroid-only
    [centroid (struct->vector centroid)]
    [centroid map-structures])
  ))

(defn fold
  "Folds a sequence of RNA and returns only the target
   structure. Target structure can either be centroid or MFE."
  
  [s & {:keys [foldtype]
        :or {foldtype "mfe"}}]
  (cond
   (= foldtype "mfe")
   (->> ((shell/sh "RNAfold"
                   "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                   "--noPS"
                   :in s)
         :out)
        (str/split-lines)
        second
        (str/split #" ")
        first)
   
   (= foldtype "centroid")
   (first (suboptimals s 10000))))

(defn neutrality
  "takes a string s and  returns neutrality of the seq  when compared
   to each of the 3L 1-mutant neighbors"

  [s & {foldtype :foldtype cons :cons n :n
        :or {foldtype "mfe" cons nil
             n (+ 3 (count
                     (filter #(> % 0.05)
                             (map (fn [i]
                                    (CREl i (fold s) :alpha ["(" "." ")"] :limit 20))
                                  (range 4 20)))))}}]
  (let [st (fold s foldtype)
        neighbors (mutant-neighbor s)
        neutfn (fn [x] (/ (- (count s) x)
                       (count s)))
        ]
    (if cons
      (let [cons cons #_(->> ((profile (read-sto cons)) :structure)
                      first
                      (str/replace-re #"\:|\-" "." ))] 
        (map stats/mean
             (transpose
              (map (fn [neighbor]
                     (let [js (jensen-shannon (probs n cons) (probs n (fold neighbor)))
                           norm (levenshtein st (fold neighbor))
                           c (apply levenshtein (-> (sw cons (fold neighbor))
                                                    first
                                                    rest))
                           rnadist (->> ((shell/sh "RNAdistance" 
                                                   :in (str cons "\n" (fold neighbor)))
                                         :out)
                                        (re-find #"\d+" )
                                        (Integer/parseInt))
                           ]
                       [(neutfn norm) (neutfn c) (neutfn rnadist) (- 1 js)]))
                   neighbors))))
      (stats/mean
       (map (fn [neighbor]
              (neutfn (levenshtein st (fold neighbor))))
            neighbors))
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


(defn subopt-overlap-seq
  "Determine the percent overlap of each suboptimal structure of a
  sequence s to the consensus structure. compares against n suboptimal
  structures.
  returns a map of frequencies where k=%overlap and v=frequency."

  [s cons-keys n]
  (let [[cent substruct] (suboptimals s n :centroid-only false)]
                          ;;takes percent overlap and
                          ;;reduces it to a freqmap to
                          ;;save memeory
    (frequencies (map (fn [ks]
                        ;;percent overlap
                        (/ (count (sets/intersection cons-keys
                                                     (set (keys ks))))
                           (count cons-keys)))
                      substruct))))

(defn subopt-overlap-sto
  "Takes a sto file and finds the suboptimal structure for the WT and
   each of its 1-mutant neighbors. Returns a vector [sto
   list-of-lists-of-maps] where each list-of-maps is the percent
   overlap for a particular sequence and its 1-mutant neighbors. The
   maps are the particular percent overlap for a 1-mutant neighbor."
  
  [sto & {:keys [ncore nsubopt altname]
          :or {ncore 6 nsubopt 1000 altname sto}}]
  (let [{l :seqs cons :cons} (read-sto sto :with-names true)
        cons (change-parens (first cons))]
    [altname ;return [filename data]
     (doall
      ;;go over each seq in the alignment
      (pxmap
       (fn [[nm s]] 
         (let [[s st] (remove-gaps s cons)
               neighbors (mutant-neighbor s)
               cons-keys (set (keys (struct->matrix st)))]
           ;;finds 1000 suboptimal structures and
           ;;finds the percent overlap of
           ;;suboptimal structures to the cons struct
           (doall                 
            (map (fn [neighbor]
                   ;;a freqmap of % overlap for each neighbor
                   (subopt-overlap-seq neighbor cons-keys nsubopt)) 
                 (concat (list s) neighbors)))));first element is WT rest are mut neighbors
       ncore
       l))] ;l=list of seqs in the sto
    ))

(defn avg-overlap
  "Takes the map of percent overlaps where it is organized k=file name
   and v=list of lists of frequency maps of percent overlap of 1000 suboptimal
   structures for the WT and each of its 1-mutant neighbors"
  
  [map-of-per-overlaps]
  (reduce (fn [m [k list-lists-maps]]
            (let [list-map (->> list-lists-maps 
                                (apply concat) ;combines data from all sequences
                                (apply merge-with +)) ;combines the
                                        ;freqmaps into 1 map
                  avg (double (mean list-map))
                  sd (double (sd list-map))
                  med (double (median-est list-map))]
              (assoc m k [med avg sd])))
          {} map-of-per-overlaps))

(defn main
  "Main function for determining neutrality of all sequences and their
  1-mutant neighbors by finding the percent overlap between the
  consensus structure and its suboptimal structures.
  :pos is a boolean."

  [& {:keys [pos]
      :or {pos true}}]
  (let [fdir (str homedir "/bin/gaisr/trainset2/"
                  (if :pos "pos/" "neg/"))
        fs (str/split-lines ((shell/sh "ls" :dir fdir) :out))
        fsto (filter #(re-find #"\.sto" %) fs)]
    (doall
     ;;go over each sto in the dir provided
     (for [sto (take 3 (filter #(re-seq #"RF00555-seed" %) fsto))]
       (subopt-overlap-sto sto)
       ))))


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

#_(defn evaluate-inverse-cons4
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


(comment 
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


(timefn (fn [] 
                      (let [fdir "/home/peis/bin/gaisr/trainset2/neg/"
                            fs (str/split-lines ((shell/sh "ls" :dir fdir) :out))
                            fsto (filter #(re-find #"\.sto" %) fs)
                            ffasta (map #(str (str/butlast 3 %) "fasta") fsto)
                            ]
                        (doall 
                         (map (fn [sto fasta]
                                (let [l (->> (io/read-lines (str fdir fasta))
                                             rest
                                             (take-nth 2)
                                             )]
                                  [sto
                                   (doall (pxmap
                                    (fn [s] (neutrality s :foldtype "centroid" :cons (str fdir sto)))
                                    6
                                    l))]
                                  ))
                              fsto ffasta)))))


(let [negs (into {} (first foo))
                  pos (into {} neut-pos)
                  get-vals (fn [m] (filter #(re-find #"RF00559" %) (keys m)))
                  x (flatten (map (fn [k]
                                    (map last (get negs k)))
                                  (get-vals negs)))]
              (map (fn [wt]
                     [(count (remove #(< % wt) x)) (count x)])
                   (map last (pos (first (get-vals pos))))))

;;;calculate neutrality of all stos in a folder. uses centroid
;;;structure and compares to consensus structure after adjusted to the
;;;length of the sequence s.
(def foo (future (timefn (fn [] 
                      (let [fdir "/home/peis/bin/gaisr/trainset2/pos/"
                            fs (str/split-lines ((shell/sh "ls" :dir fdir) :out))
                            fsto (filter #(re-find #"\.sto" %) fs)
                            
                            ]
                        (doall 
                         (map (fn [sto]
                                  (let [{l :seqs cons :cons} (read-sto (str fdir sto) :with-names true)
                                        cons (change-parens (first cons))]
                                  [sto
                                   (doall (pxmap
                                           (fn [[nm s]] (let [[s st] (remove-gaps s cons)] 
                                                     (neutrality s :foldtype "centroid" :cons st)))
                                    6
                                    l))]
                                  ))
                              (drop 1 fsto))))))))

;;;figure out the loop/bulge/stem lengths averaged across RFAM
(let [fdir "/data2/Bio/RFAM/"
      fs (filter #(re-find #"\.sto" %) 
                 (str/split-lines ((shell/sh "ls" :dir fdir) :out)))
      m (map (fn [sto]
               (let [f (read-sto (str fdir sto))
                     cons (change-parens (first (f :cons)))
                     loop? (fn [st]
                             (re-seq #"\([\.AaBb]+?\)" st)) ;loops have open close
                     bulge? (fn [st]
                              (concat (re-seq #"\([\.AaBb]+?\(" st)  ;find open paren buldges
                                      (re-seq #"\)[\.AaBb]+?\)" st)));find close paren buldges
                     ;;finds the stem as long as there is some
                     ;;pattern of 0-2 gap chars in the stem
                     ;;surrounded by open or close parens
                     stem? (fn [st]
                             (concat (re-seq #"\({2,}+" st)   
                                                 (re-seq #"\){2,}+" st))
                             #_(concat (re-seq #"\(+[\.AaBb]{0,2}\(*[\.AaBb]{0,2}\(+" st)   
                                     (re-seq #"\)+[\.AaBb]{0,2}\)*[\.AaBb]{0,2}\)+" st)))
                     ;;stems count all parts of the regex
                     ;;otherwise the open/close parens used in
                     ;;the regex are not counted thus -2 from
                     ;;the count
                     length? (fn [col & {stem :stem
                                        :or {stem false}}]
                               (if stem
                                 (map count col)
                                 (map #(- (count %) 2) col)))]
                 [(length? (loop? cons))
                  (length? (bulge? cons))
                  (length? (stem? cons) :stem true)]))
             fs)]
  (map #(vector (stats/mean (apply concat %))
         (stats/sd (apply concat %)))
       (transpose m
                  )))

;;;testing/looking at structure difference. this has revealed that
;;;levenshtein distance is not a very good way to compare the
;;;difference. maybe switch to RNAdistance since the structures are
;;;being shortened according to the sequence length.
(timefn (fn [] 
                      (let [fdir "/home/peis/bin/gaisr/trainset2/pos/"
                            fs (str/split-lines ((shell/sh "ls" :dir fdir) :out))
                            fsto (filter #(re-find #"\.sto" %) fs)
                            
                            ]
                        (doall 
                         (map (fn [sto]
                                  (let [{l :seqs cons :cons} (read-sto (str fdir sto) :with-names true)
                                        cons (change-parens (first cons))]
                                  [sto
                                   (doall (map
                                           (fn [[nm s]] (let [[s st] (remove-gaps s cons)
                                                             centroid (first (suboptimals s 10000)) 
                                                             n 9 l (count s)]
                                                         (prn "cent" centroid)
                                                         (prn "cons" st)
                                                         (prn (- 1 (jensen-shannon (probs n cons) (probs n centroid))) 
                                                              (- 1 (double (/ (levenshtein cons centroid) l)))
                                                              (- 1 (double (/ (->> ((shell/sh "RNAdistance" 
                                                                                         :in (str cons "\n" centroid))
                                                                               :out)
                                                                              (re-find #"\d+" )
                                                                              (Integer/parseInt)) 
                                                                         l))))))
                                    
                                    l))]
                                  ))
                              ["RF00555-seed.4.sto" "RF00558-seed.3.sto"
                               "RF00559-seed.7.sto" "RF00167-seed.9.sto"])))))





;;;plotting the data from abc
(let [lines (map (fn [ms] 
                   (let [x (map (fn [m] (freqn->list m)) ms) ;turn each
                                        ;1-mut neighbor map into list
                         wt (stats/mean (first x)) ;first one is WT
                         ys (map #(stats/mean (flatten %)) (partition-all 3 (rest x)))] ;groups remaining lists into groups of 3 (3 1-mut neighbor per position) and averages all values
                     (cons wt ys))) 
                 (get abcdneg "RF00555-seed.neg1.4.sto"))
      l (charts/xy-plot (range 200) (first lines) :title "neg RF00555.4" :series 1 :legend true
                        :x-label "position" :y-label "mut % overlap with cons")]
  (view l) (map (fn [i y]
                  (charts/add-lines l (range 200) y :series-label i))
                (iterate inc 2) (rest lines)))

;;;stuff to process the resulting def above 'abc'
(io/with-out-writer "/home/peis/bin/gaisr/src/temp.txt"
  (doseq [i (map (fn [x] (map double x))
                 (->> abc first first second first ))]
    (println (str/join "," i))))

;;;R code to process the pvals of the data

;  > data<-read.csv("~/bin/gaisr/src/temp.txt",header=FALSE)
;  > result<-apply(data, 1, function(y) { ht<-t.test(data[1, ],y,conf.level=0.99); ht$p.value})
;  > write.csv(result,"~/bin/gaisr/src/out-pval.csv")


;;;read in resulting pvals from R output
(filter #(<= (second %) 0.01)
        (map (fn [i]
               (vec (map read-string (->> i (str/split #",")))))
             (drop 1 (io/read-lines "src/out-pval.csv"))))
                                      
;;;look at data from abc (abcd is a map made of abc). groups all of
;;;hte data together for each sto to make a histogram to see what the
;;;overall distribution of the data is.
(view (charts/histogram (flatten (map (fn[[x n]]
                                        (repeat n x))(->> (get abcd "RF00555-seed.4.sto") (map #(apply merge-with + %)) (apply merge-with +))))))


;;;make a graph of the avg overlap of all 1-mut-structures (3) with
;;;cons-structure at a position i from output abc
(let [x (map (fn [m] (freqn->list m)) (->> foo first first second second))
                  wt (stats/mean (first x))
                  ys (map #(stats/mean (flatten %)) (partition-all 3 (rest x)))]
  (view (charts/line-chart (range (count x)) (cons wt ys)
                           :title "neg RF00555"
                           :x-label "position"
                           :y-label "percent overlap with cons structure")))

;;;get the average percent overlap for each sto in the training set
(def foo (let [pos (for [k (keys abcd)] ;loop over all keys
                     ;;the average comes from all %overlap values in
                     ;;the sto; therefore all points have equal weight
                          (stats/mean (flatten (map (fn [ms] 
                                 (let [x (map (fn [m] (freqn->list m)) ms) ;convert a freq map to a list
                                       ] 
                                   x))
                               (abcd k)))))
                  neg (for [k (keys abcdneg)]
                          (stats/mean (flatten (map (fn [ms] 
                                 (let [x (map (fn [m] (freqn->list m)) ms)
                                       ] 
                                   x))
                               (abcdneg k)))))]
           [pos neg]))


(defn avg-overlap-old 
  "DEPRECATED
   takes the map of percent overlaps where it is organized k=file name
   and v=list of lists of maps of percent overlap of 1000 suboptimal
   structures for the WT and each of its 1-mutant neighbors"

  [map-of-per-overlaps]
  (reduce (fn [m [k list-lists-maps]]
            (let [sum-freqmap (fn ;;figures out the sum of a freqmap
                                ;;where the k=value and v=number of appearance
                                [m]
                                (reduce (fn [x [val n]]
                                          [(+ (* n val) (first x)) ;sum (v * n)
                                           (+ n (second x))])
                                        [0 0] m))
                  ;;finds the sums for each sequence in the alignment
                  sums (map (fn [list-maps] 
                              (let [x (transpose (map sum-freqmap list-maps)) ;convert a freq map a vector [total-of-values total-number-counted]
                                    ]
                                ;;transpose makes all total-of-values appear together
                                [(sum (first x)) 
                                 (sum (second x))]))
                            list-lists-maps)
                  ;;adds the sums for percent overlap  found for each sequence in the alignment
                  avg (let [totals (map sum (transpose sums))] ;return average
                        (double (/ (first totals) (second totals))))
                  sd (let [sqs (map (fn [[val n]]
                                     [(* (- val avg)
                                         (- val avg))
                                      n])
                                   (transpose sums))]
                           (-> (sum (first sqs))
                               (/ (sum (second sqs)))
                               Math/sqrt)) ]
              (assoc m k [avg sd])))
          {} map-of-per-overlaps))

(defn overlap-per-seq
  "Takes an entry from teh map-of-lists-of-lists-of-maps data
   structure. This entry has a format of [nm data]. It averages the
   %overlap for all mutations at a position. Returns a vector of
   average %overlap by position."

  [per-overlaps]
  (let [[nm data] per-overlaps]
    (map (fn [ea-seq]
           (let [wt (mean (first ea-seq))
                 ;;merges the 3 mutations at a position and finds the mean
                 muts (for [mut (partition-all 3 (rest ea-seq))]
                        (->> mut
                             (apply merge-with +) ;merge freqmaps together
                             mean))]
             (vec (cons wt muts)))) ;returns vector starting with wt
                                    ;then muts
         data)))

(defn chart-overlap-sto
  "Takes an entry from the map-of-lists-of-lists-of-maps data
   structure. It calls the overlap-per-seq to find the avarge
   %overlaps at each position. Returns a graph of %overlap for each
   sto. Each line on graph represents 1 sequence from the sto."

  [per-overlaps]
  (let [lines (overlap-per-seq per-overlaps)
        l (charts/xy-plot (range 200) (first lines) :title "neg RF00555.4" :series 1 :legend true
                          :x-label "position" :y-label "mut % overlap with cons")]
    (view l)
    (map (fn [i y]
           (charts/add-lines l (range 200) y :series-label i))
         (iterate inc 2) (rest lines))))

(defn subopt-significance
  "Takes an input sto and estimates the significance of the suboptimal
   structure overlap when compared to n random simulated
   alignments. Returns the map-of-lists-of-lists-of-maps."
  
  [insto & {:keys [ncores nsamples]
            :or {ncores 3 nsamples 3}}]
  (concat [(subopt-overlap-sto insto :altname (fs/basename insto))] ;wt subopt overlap
          ;;shuffled version of sto
          (pxmap (fn [randsto]                   
                   (let [result (subopt-overlap-sto randsto :altname (fs/basename randsto))] ;find subopt overlap
                     (fs/rm randsto)
                     result))
                 ncores
                 (repeatedly nsamples #(sto->randsto insto (fs/tempfile))) ;create n stos
                 )))

(def foo (let [insto "/home/peis/bin/gaisr/trainset2/pos/RF00167-seed.4.sto"] 
           (subopt-significance insto)))


)


