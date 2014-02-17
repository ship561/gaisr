(ns robustness.sandbox
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            ;;[incanter.stats :as stats]
            [incanter.charts :as charts]
            ;;[clojure.contrib.json :as json]
            [clojure.set :as sets]
            [clojure.core.reducers :as r]
            [edu.bc.utils.snippets-utils]
            [edu.bc.bio.sequtils.alignment-info])
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        [edu.bc.utils.snippets-math :only (median-est)]
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens sto->randsto read-clj
                         parse-dotps)]
        edu.bc.utils.fold-ops
        [incanter.core :only (view)]
        invfold
        refold
        smith-waterman
        ))

(timefn (fn [] 
          (let [fdir "/home/peis/bin/gaisr/trainset2/"
                fs (str/split-lines ((shell/sh "ls" :dir (str fdir "pos/")) :out))
                fsto fs
                ffasta (map #(str (str/butlast 3 %) "fasta") fs)
                x (fn [s cons n] (let [cons (->> ((read-sto cons) :cons)
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
                                  [(pearson-correlation (first tx) (second tx))
                                   (pearson-correlation (first tx) (third tx))
                                   (pearson-correlation (second tx) (third tx))]))
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
                                     (let [{l :seqs cons :cons} (read-sto (str fdir sto) :info :both)
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
                    (let [{l :seqs cons :cons} (read-sto (str fdir sto) :info :both)
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





(def foo (let [insto "/home/peis/bin/gaisr/trainset2/pos/RF00167-seed.4.sto"] 
           (subopt-significance insto)))

(def signif (future (timefn (fn [] (let [fdir "/home/peis/bin/gaisr/trainset2/pos/"]
                                    (doall (map (fn [insto]
                                                  [insto (into {} (subopt-significance (str fdir insto) :ncores 2 :nsamples 100))])
                                                (filter #(re-find #"\.3\.sto" %) (fs/listdir fdir)))))))))

;;;show robustness of structures when comparing a sequence with its cons structure and the
;;;inverse folded sequence which is similar to the cons structure. This should show that the neutrality
;;;of the original sequence is significantly robust when compared to its inverse folded sequences.
(let [sto "/home/kitia/bin/gaisr/trainset2/pos/RF00555-seed.1.sto"
      {l :seqs cons :cons} (read-sto sto :info :both)
      cons (change-parens (first cons))]
  [sto 
   (pxmap (fn [[nm s]]
            (let [[s st] (remove-gaps s cons)
                  inv-fold (fn [target n]
                             (remove nil?
                                     (flatten 
                                      (map (fn [[s ensemble]]
                                             (re-find #"\w+" s))
                                           (->> ((shell/sh "RNAinverse"
                                                           "-Fmp"
                                                           (str "-R" n)
                                                           "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                                           :in target)
                                                 :out)
                                                (str/split-lines)
                                                (partition-all 2))
                                           ))))
                  inv-seq (inv-fold st 100)
                  neut (map (fn [x]
                              (subopt-overlap-neighbors x st :nsubopt 1000))
                            (concat (list s) inv-seq))]
              (map (fn [x] (->> x (apply merge-with +) mean double)) neut)))
          2
          (take 1 l))])

;;;analyze resulting data for avg neutrality
(let [[wt & muts]  (->> foo second transpose )]
  [(stats/mean wt) (stats/mean (flatten muts))])


(let [s "UCAAACAAUGAGAACAUUACUUAUUUAUGUCACGAA...UGGGCGUGACGUUUCUACAAGGUG.CCGU.AA.CACCUAACAAUAAGUAAGCUAAUUUAGUCA"
      st "...............((((((((((...(.((((.........)))).)(((......((((((.......)))))))))))))))).....)))......."
      [s st] (remove-gaps s st)
      cons-keys (set (keys (struct->matrix st)))
      n 1000
      neighbors (mutant-neighbor s :with-names true)
      wt (map #(probs 1 %) (transpose (fold s :foldtype "RNAsubopt" :n n)))
      loi (->> (for [[nm neighbor] neighbors
                     :let [mut (map #(probs 1 %) (transpose (fold neighbor :foldtype "RNAsubopt" :n n)))
                           overlap (double (mean (subopt-overlap-seq neighbor cons-keys n)))]]
                 (map (fn [c1 c2]
                        [nm (jensen-shannon c1 c2) overlap])
                      wt mut))
               transpose
               (map (fn [jsd]
                      (remove #(or (< (second %) 0.9)
                                   (> (third %) 0.3))
                              jsd))
                    )
               (interleave (iterate inc 0) )
               (partition-all 2 )
               (remove #(empty? (second %)) )
               (map vec ))
      neighbors (into {} neighbors)
      ]
  (io/with-out-writer "/home/kitia/bin/gaisr/robustness/struct-confirm.txt"
    (doseq [[pos x] loi
            [mutnm h overlap] x]
      (prn "wt" pos)
      (prn s)
      (prn (apply str (repeat 10 "0123456789")))
      (doseq [i (fold s :foldtype "RNAsubopt" :n 3)]
        (prn i))
      (prn mutnm pos h overlap)
      (prn (neighbors mutnm))
      (prn (apply str (repeat 10 "0123456789")))
      (doseq [i (fold (neighbors mutnm) :foldtype "RNAsubopt" :n 3)]
        (prn i)))))


(let [sto "/home/kitia/bin/gaisr/trainset2/pos/RF00555-seed.1.sto"
      {l :seqs cons :cons} (read-sto sto :info :both)
      cons (change-parens (first cons))]
  (io/with-out-writer "/home/kitia/bin/gaisr/robustness/struct-confirm.txt"
    (prn sto)
    (doseq [[seqnm s] l]
      (let [[s st] (remove-gaps s cons)
            cons-keys (set (keys (struct->matrix st)))
            n 1000
            neighbors (mutant-neighbor s :with-names true)
            wt (map #(probs 1 %) (transpose (fold s :foldtype "RNAsubopt" :n n)))
            loi (->> (for [[nm neighbor] neighbors
                           :let [mut (map #(probs 1 %) (transpose (fold neighbor :foldtype "RNAsubopt" :n n)))
                                 overlap (double (mean (subopt-overlap-seq neighbor cons-keys n)))]]
                       (map (fn [c1 c2]
                              [nm (jensen-shannon c1 c2) overlap])
                            wt mut))
                     transpose
                     (map (fn [jsd]
                            (remove #(or (< (second %) 0.8)
                                         (> (third %) 0.4))
                                    jsd))
                          )
                     (interleave (iterate inc 0) )
                     (partition-all 2 )
                     (remove #(empty? (second %)) )
                     (map vec ))
            neighbors (into {} neighbors)
            ]
        (prn seqnm)
        (doseq [[pos x] loi
                [mutnm h overlap] x]
          (prn "wt" pos)
          (prn s)
          (prn (apply str (repeat 10 "0123456789")))
          (doseq [i (fold s :foldtype "RNAsubopt" :n 3)]
            (prn i))
          (prn mutnm pos h overlap)
          (prn (neighbors mutnm))
          (prn (apply str (repeat 10 "0123456789")))
          (doseq [i (fold (neighbors mutnm) :foldtype "RNAsubopt" :n 3)]
            (prn i)))))))


(let [s "UCAAACAAUGAGAACAUUACUUAUUUAUGUCACGAA...UGGGCGUGACGUUUCUACAAGGUG.CCGU.AA.CACCUAACAAUAAGUAAGCUAAUUUAGUCA"
      st "...............((((((((((...(.((((.........)))).)(((......((((((.......)))))))))))))))).....)))......."
      [s st] (remove-gaps s st)
      cons-keys (set (keys (struct->matrix st)))
      n 1000
      neighbors (into {} (mutant-neighbor s :with-names true))
      wt (map #(probs 1 %) (transpose (fold s :foldtype "RNAsubopt" :n n)))
      loi (->> (for [[nm neighbor] neighbors
                     :let [mut (map #(probs 1 %) (transpose (fold neighbor :foldtype "RNAsubopt" :n n)))
                           overlap (double (mean (subopt-overlap-seq neighbor cons-keys n)))]]
                 (map (fn [i c1 c2]
                        [nm [i (jensen-shannon c1 c2) overlap]])
                      (iterate inc 0) wt mut))
               (map (fn [jsd]
                      (remove #(or (< (-> % second second) 0.9)
                                   (> (-> % second third) 0.3))
                              jsd))
                    )
               (remove #(empty? %) )
               (map vec ))
      ]
  (io/with-out-writer "/home/kitia/bin/gaisr/robustness/struct-confirm.txt"
    (doseq [x loi
            [mutnm [pos h overlap]] x]
      (prn "wt" pos)
      (prn s)
      (prn (apply str (repeat 10 "0123456789")))
      (doseq [i (fold s :foldtype "RNAsubopt" :n 3)]
        (prn i))
      (prn mutnm pos h overlap)
      (prn (neighbors mutnm))
      (prn (apply str (repeat 10 "0123456789")))
      (doseq [i (fold (neighbors mutnm) :foldtype "RNAsubopt" :n 3)]
        (prn i)))))

(let [sto "/home/kitia/bin/gaisr/trainset2/pos/RF00555-seed.1.sto"
      {l :seqs cons :cons} (read-sto sto :info :both)
      cons (change-parens (first cons))]
  (io/with-out-writer "/home/kitia/bin/gaisr/robustness/struct-confirm.txt"
    (prn sto)
    (doseq [[seqnm s] l]
      (let [[s st] (remove-gaps s cons)
            cons-keys (set (keys (struct->matrix st)))
            n 1000
            neighbors (mutant-neighbor s :with-names true)
            wt (map #(probs 1 %) (transpose (fold s :foldtype "RNAsubopt" :n n)))
            loi (->> (for [[nm neighbor] neighbors
                           :let [mut (map #(probs 1 %) (transpose (fold neighbor :foldtype "RNAsubopt" :n n)))
                                 overlap (double (mean (subopt-overlap-seq neighbor cons-keys n)))]]
                       (map (fn [i c1 c2]
                              [nm [i (jensen-shannon c1 c2) overlap]])
                            (iterate inc 0) wt mut))
                     identity;transpose
                     (map (fn [jsd]
                            (remove #(or (< (-> % second second) 0.8)
                                         (> (-> % second third) 0.4))
                                    jsd))
                          )
                     (remove #(empty? %) )
                     (map vec ))
            neighbors (into {} neighbors)
            ]
        (prn seqnm)
        (doseq [x loi
                [mutnm [pos h overlap]] x]
          (prn "wt" pos)
          (prn s)
          (prn (apply str (repeat 10 "0123456789")))
          (doseq [i (fold s :foldtype "RNAsubopt" :n 3)]
            (prn i))
          (prn mutnm pos h overlap)
          (prn (neighbors mutnm))
          (prn (apply str (repeat 10 "0123456789")))
          (doseq [i (fold (neighbors mutnm) :foldtype "RNAsubopt" :n 3)]
            (prn i)))))))


(let [fdir "/home/kitia/bin/gaisr/trainset2/pos/"
      invseqs (filter #(re-find #"\.inv.clj" %) (fs/listdir fdir))
      stos (map #(str (str/butlast 7 %) "sto") invseqs)
      foo (map (fn [sto invseq]
                 (let [sto-seq (->> (read-sto (str fdir sto) :info :both) :seqs (into {}))
                       inv-seqs (->> (read-clj (str fdir invseq)) (into {}))]
                   (vec (cons sto 
                              (for [[nm s] sto-seq] 
                                (cons (GC-content (str/replace-re #"\." "" s))
                                      (map GC-content (get inv-seqs nm))))))))
               stos invseqs)]
  (prn (count foo))(->> (drop 11 foo)
                        first
                        rest
                        (map #(view (charts/histogram %)) )))





;;;for finding the bpdistance
(let [fun-sto (fn [sto distfun]
                (let [{l :seqs cons :cons} (read-sto sto :info :both)
                      cons (change-parens (first cons))
                      altname sto]
                  [altname ;return [filename data]
                   (doall
                    ;;go over each seq in the alignment
                    (map
                     (fn [[nm s]] 
                       (let [[s st cons-keys] (degap-conskeys s cons)]
                         ;;finds 1000 suboptimal structures and
                         ;;finds the percent overlap of
                         ;;suboptimal structures to the cons struct
                         (distfun s st)))
                                        ;ncore
                     l))] ;l=list of seqs in the sto
                  ))
      main (fn [& args] 
             (let [[opts _ usage] (cli args
                                       ["-f" "--file" "file(s) to check neutrality for"
                                        :parse-fn #(str/split #" " %) ;create list of files
                                        :default nil]
                                       ["-D" "--distfn" "distance function to use" :default bpsomething]
                                       ["-o" "--outfile" "file to write to" :default nil]
                                       ["-di" "--dir" "dir in which files are located" 
                                        :default nil]
                                       ["-nc" "--ncore" "number cores to use" 
                                        :parse-fn #(Integer/parseInt %) :default 6]
                                       ["-h" "--help" "usage" :default nil :flag true])
                   stos (if (opts :dir)
                          (map #(fs/join (opts :dir) %) (opts :file))
                          (opts :file))
                   neutrality (pxmap (fn [sto] ;loop over stos
                                       (prn :workingon (fs/basename sto))
                                       #_(fun-sto sto (opts :distfn))
                                       (generic-dist-sto sto (opts :distfn)))
                                     20
                                     stos)]
               (cond
                (or (nil? args) (opts :help)) (print usage) ;usage help
                (not (nil? (opts :outfile))) (io/with-out-writer (opts :outfile)
                                               (prn (vec neutrality))) ;data to file
                :else
                (doall neutrality))))
      ]
  (main "-f" (->> (fs/directory-files "/home/peis/bin/gaisr/trainset3/neg/5prime/" ".sto")
                  (filter #(> (-> % edu.bc.bio.sequtils.alignment-info/alignment-quality first) 0.05) )
                  (str/join " " ))
        "-D" #(bpsomething %1 %2 :bp true)))

(edu.bc.utils.snippets-utils/with-out-appender
  "/home/peis/bin/gaisr/robustness/compare-neutrality/neutrality-distribution-trainset3.csv" 
  (doseq [[fnm vs] @bpdist-3prime-neg2
          v vs
          :let [sto (fs/basename fnm)
                {nm :name type :type} (edu.bc.bio.sequtils.alignment-info/parse-sto-function (str (re-find #"RF\d+\-seed" sto) ".sto"))]]
    (println (str sto "," v ",bpdist-3prime-neg2," nm "," (apply str type)))))

#_(main "-f" (->> (fs/directory-files "/home/peis/bin/gaisr/trainset3/neg/3prime/" ".sto")
                  (filter valid-seq-struct )
                  (map fs/basename)
                  (str/join " " ))
        "--dir" "/home/peis/bin/gaisr/trainset3/neg/3prime/"
        "-D" #(pccsomething %1 %2))


(io/with-out-writer "/home/peis/bin/gaisr/robustness/compare-bpdist12.clj"
  (println ";;;compares the WT structure to the centroid structure of the 1-mut neighbors. turns the structures into a vector of 0's and 1's where 0 is unpaired 1 is paired. Uses the pearsons correlation coefficient to find the similarity between the 2 structures. the equation 1-0.5*(1-pcc) is used to normalize the data to [0 1]. uses pccsomething over the 5' NEGative training set.")
  (prn (mapv #(vector (fs/basename (first %))
                      (vec (second %)))
             @bar)))

;;random stuff to add groupings to the csv. temp file created. must
;;delete afterwards. 
(let [[header & data] (csv/parse-csv (slurp "/home/peis/bin/gaisr/robustness/compare-neutrality/neutrality-distribution-trainset3.csv"))
      update-header (conj header "distance metric")]
  (io/with-out-writer "/home/peis/bin/gaisr/robustness/compare-neutrality/neutrality-distribution-trainset3-1.csv"
    (print (csv/write-csv (conj (map (fn [d] 
                                       (conj d (re-find #"^\w+" (third d))))
                                     data)
                                update-header)))))

(defn data->csv
  "takes the data from main-subopt-overlap and makes it into csv
  format and appends it to the outcsv file. example

  data = @bpdist-3prime-neg2
  outcsv = \"/home/peis/bin/gaisr/robustness/compare-neutrality/neutrality-distribution-trainset3.csv\"
  dist-measure =  \"bpdist-3prime-neg2\" "
  
  [data outcsv dist-measure]
  (let [outfn (fn []
                (doseq [[fnm vs] data
                        v (map str vs)
                        :let [sto (fs/basename fnm)
                              {nm :name type :type} (edu.bc.bio.sequtils.alignment-info/parse-sto-function (str (re-find #"RF\d+\-seed" sto) ".sto"))
                              dist-metric (re-find #"^\w+" dist-measure)]]
                  (print (csv/write-csv [[sto v dist-measure nm (apply str type) dist-metric]]))))]
    (if (fs/exists? outcsv)
      (edu.bc.utils.snippets-utils/with-out-appender outcsv (outfn))
      (io/with-out-writer outcsv
        (println "file name,neutrality,distance measure,ID,function,distance metric")
        (outfn)))))







;;;finding clusters of structures
(time (let [n 1000
                          s "AUAGCGUACAAAAAAAUCCCUCAUAUAAUUUUGGGAAUAUGGCCCAAAAGUCUCUACCCAAUAACCGUAAAUUAUUGGACUAUGAGGGAAAGGAUCUGUUUU"
                          sts (as-> (fold s {:foldmethod :RNAsubopt :n n}) sts 
                                    (interleave  sts (range n))
                                    (partition 2 sts) 
                                    (vec sts))
                          dm (r/fold 50 (fn ([] {})
                                       ([l r] (merge l r)))
                                     (fn ([] {})
                                       ([M [sti i]]
                                          (as-> sts stj
                                                (r/filter (fn [[_ j]] (< i j)) stj)
                                                (r/map (fn [[stj j]]
                                                         [[i j](bpdist sti stj :bpdist true)])
                                                       stj)
                                                (into {} stj)
                                                (merge M stj))))
                                     sts)]
                            
                      (io/with-out-writer "/home/kitia/TMP/dist.matrix"
                        (doseq [i (range n)]
                          (->> (map (fn [j]
                                      (if (= i j)
                                        0
                                        (get dm [i j] (dm [j i])))) ;dist
                                    (range n))
                               (str/join "," )
                               print)
                          (println "")))))
