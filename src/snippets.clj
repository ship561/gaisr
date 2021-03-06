(ns snippets
  (:require [clojure.contrib.string :as str]
            [clojure.java.shell :as shell]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs]
            [incanter.stats :as stats]
            [incanter.charts :as charts]
            [clojure.set :as sets])
  (:use edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        ;smith_waterman
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.core :only (view)]
        refold))

(def f
  (let [info (mutual_info (profile (read-sto "/home/kitia/Downloads/S15_101711UBedit.sto")))
        x (reduce (fn [m [k v]]
                    (assoc m k v))
                  {} info)]
    x))

(view (heat-map #(if (and (integer? %1) (integer? %2))
                   (get f [%1 %2] 0)
                   0)
                0 200 0 200))

(defn randseq
  "takes a sequence and toggles changes so that it will become a
   random sequence that shares a thr of sequence identity"

  [inseq thr]
  (let [b #{"A" "C" "G" "U"}
        change? (fn [base]
                  (if (> (rand) thr)
                    (get (vec (disj b base)) (rand-int 3))
                    base))
        energy (fn [s]
                 (Double/parseDouble
                  (re-find #"\-*\d*.\d+"
                           (-> (str s "\n")
                               ((fn [x] (shell/sh "echo" x)))
                               ((fn [x] (shell/sh "RNAfold" "--noPS" "-P" "/home/kitia/Desktop/ViennaRNA-2.0.1/rna_andronescu2007.par" :in x)))))))]
    (energy (apply str
           (for [i (rest (str/split #"" inseq))]
             (change? i))))))


(defn markov-chain-seq
  "Uses a MCMC to make a new sequence of similar dinucleotide
   frequency. Must give an input sequence s. The function will find
   the conditional probability using P(X|Y)= P(X,Y)/P(Y). If a
   dinucleotide is not in the sequence, then 1e-6 is used as a
   pseudocount to give a small probabilty of getting it."

  [s]
  (let [example "AAGGUACCUU"
        s (str/replace-re #"T" "U" s)
        mono (reduce (fn [m [k v]]
                       (assoc m (str k) v)) 
                     {} (probs 1 s))
        di (probs 2 s)
        conditional (fn [cur-base]
                      (reduce (fn [m y]
                                (assoc m y (/ (get di (str cur-base y) 0.000001)
                                              (get mono cur-base 0.000001))))
                              {} #{"A" "C" "G" "U"}))
        normalize (fn [m]
                    (reduce (fn [a [k v]]
                              (assoc a k (/ v (sum m))))
                            {} m))
        ]
    (loop [i (range 1 (count s))
           ss (str (first s))]
      ;;(prn (first i) ss)
      (if-not (empty? i)
        (let [r (rand)
              markovstep (into {}
                               (map (fn [k v]
                                      [k v])
                                    (keys (normalize (conditional (str (last ss)))))
                                    (drop 1 (reductions (fn [c v] 
                                                          (+ c v))
                                                        0 (vals (normalize (conditional (str (last ss)))))))))]
          ;;(prn r markovstep) 
          (recur (rest i)
                 (cond 
                  (< r (-> markovstep first val))
                  (str ss (-> markovstep first key))
                  (< r (-> markovstep second val))
                  (str ss (-> markovstep second key))
                  (< r (->> markovstep (drop 2) first val))
                  (str ss (->> markovstep (drop 2) first key))
                  :else
                  (str ss (-> markovstep last key)))))
        ss))))

(defn nussinov
  "Uses the Nussinov algorithm to compute an optimal RNA structure by
  maximizing base pairs in the structure. The function requires an
  input string s. The output is a list of base pair locations [i
  j]. It will also print out the sequence and the structure so that it
  can be visually inspected. An example sequence of 'GGGAAAUCC' will
  give the answer ([2 6] [1 7] [0 8]). Locations are 0 based (ie seq
  goes from 0 to n-1)."
  
  [s]
  (let [s (.toUpperCase s)
        n (count s)
        delta (fn [i j] ;;base pairing score
                (let [bp {"AU" 1 "UA" 1 "GC" 1 "CG" 1 "GU" 1 "UG" 1}
                      b1 (subs s i (inc i))
                      b2 (subs s j (inc j))]
                  (get bp (str b1 b2) 0)))
        loc (filter #(and (< (second %) n) ;possible locations in the table to fill in
                          (> (- (second %) (first %)) 3)) ;;must be
                                                          ;;atleast 3 bases in a loop
                    (for [k (range (- n 1)) ;;fill in table diagonally
                          i (range n)]
                      [i (+ i k 1)]))
        gamma (reduce (fn [m [i j]] ;;table of values
                        (let [g (fn [i j] ;;determines where the
                                         ;;current gamma(i,j) score comes from
                                  (max (get m [(inc i) j] 0) ;;i unpaired
                                       (get m [i (dec j)] 0) ;;j unpaired
                                       (+ (get m [(inc i) (dec j)] 0) (delta i j)) ;;i j paired
                                       (apply max (if (empty? (range (inc i) j)) ;;bifurcation
                                                    '(-1)
                                                    (for [k (range (inc i) j)]
                                                      (+ (get m [i k] 0) (get m [(inc k) j] 0)))))))]
                          (assoc m [i j] (g i j))))
                      {} loc)
        ;;traceback starts at [0 n-1]
        traceback (loop [x (list [0 (dec n)]) ;;x is stack of positions to check
                         bp (list)]
                    (if-not (empty? x)
                      (let [[i j] (first x)
                            x (pop x)]
                        (if (< i j)
                          (let [ks (filter (fn [k]
                                             (= (+ (get gamma [i k] 0) (get gamma [(inc k) j] 0)) (get gamma [i j] 0)))
                                           (range (inc i) j))]
                                        ;(prn i j (get gamma [i j] 0) ks)
                            (cond
                             (= (get gamma [(inc i) j] 0) (get gamma [i j] 0)) ;;i unpaired
                             (recur (conj x [(inc i) j])
                                    bp)
                             (= (get gamma [i (dec j)] 0) (get gamma [i j] 0)) ;;j unpaired
                             (recur (conj x [i (dec j)])
                                    bp)
                             (= (+ (get gamma [(inc i) (dec j)] 0) (delta i j)) (get gamma [i j] 0)) ;;i j base paire
                             (recur (conj x [(inc i) (dec j)])
                                    (conj bp [i j]))
                             (not (empty? ks)) ;;bifurcation
                             (recur (conj x [i (first ks)] [(inc (first ks)) j])
                                    bp)
                                        ;:else ;;not sure if necessary. only 4 possible conditions
                                        ;(recur x
                                        ;       bp)
                             ))
                          (recur x 
                                 bp)))
                      bp))
        structure (loop [bp traceback
                         st (apply str (repeat n "."))]
                    (if-not (empty? bp)
                      (let [[i j] (first bp)]
                        (recur (rest bp)
                               (str (subs st 0 i) "(" (subs st (inc i) j) ")" (subs st (inc j) n))))
                      st))]
    (prn s)
    (prn structure)
    traceback))



;;code looks through list of stos and picks random sequences out of
;;the sto
(doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset/neg/list-sto-only.txt") ;;use trainset/list2 or trainset/neg/list-sto-only
        ]
  (let [c (atom 0)
        fdir "/home/kitia/bin/gaisr/trainset/neg/" ;;use either trainset/ or trainset/neg
        odir "/home/kitia/bin/gaisr/trainset2/neg/" ;;use either trainset2/pos or trainset2/neg
        [_ seq-lines _] (join-sto-fasta-lines (str fdir f) "")
        ;;f "RF01693-seed.sto"
        ]
    (if (>= (nCk (count seq-lines) 3) 10)
      (doseq [i (->> (filter #(apply distinct? %) (repeatedly #(randseqs seq-lines (+ (rand-int 4) 3))))
                     (take 60)
                     distinct
                     (take 10))]
        (prn c f)
        (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/temp.aln"
          (println "CLUSTAL W (1.83) multiple sequence alignment\n")
          (doseq [[nm sq] i]
            (cl-format true "~A~40T~A~%" nm sq)))
        (aln->sto "/home/kitia/bin/gaisr/trainset2/temp.aln" (str odir (subs f 0 (- (count f) 3)) @c ".sto"))
        (swap! c inc))
      (println f "did not work. too few sequences")
      )
    ))



(let [d "/home/kitia/bin/gaisr/trainset2/"
              ]
  (doseq [stos (fs/listdir d)
          s stos]
                     (let [fasta (str (subs s 0 (- (count s) 3)) "fasta")]
                          (sto->fasta (str d s) (str d fasta)))))



;;gets structures out of sto files and then does a variety of calcs on it
(for [f (remove #(or (= % "RF01510-seed.sto")
                     (= % "RF01510-seed.neg1.sto"))
                (io/read-lines "/home/kitia/bin/gaisr/trainset/list2.txt")) ;;use trainset/list2 or trainset/neg/list-sto-only
        ]
  (let [c (range 10)
        fdir "/home/kitia/bin/gaisr/trainset/" ;;use either trainset/ or trainset/neg
        odir "/home/kitia/bin/gaisr/trainset2/pos/" ;;use either trainset2/pos or trainset2/neg
        ;;f "RF01693-seed.sto"
        cons-lines (->> (join-sto-fasta-lines (str fdir f) "") last first second second (str/replace-re #"\<" "(" ) (str/replace-re #"\>" ")") (str/replace-re #"[^().]" "." ))
        Prcons-lines (probs 1 cons-lines)
        ]
        (catch-all (for [i c]
               (let [x (-> (join-sto-fasta-lines (str odir (subs f 0 (- (count f) 3)) i ".sto") "") last first second second)
                    [[_ s1 s2]] (sw cons-lines x)
                    s1 (str/replace-re #"\-" "." s1)
                    s2 (str/replace-re #"\-" "." s2)
                    Prx (probs 1  x)
                    Prs1 (probs 1 s1)
                    Prs2 (probs 1 s2)] #_(prn s1) #_(prn s2)
                    (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/temp.struct" (println cons-lines) (println x))
                    [(count cons-lines) (count x) (count s1) (count s2) 
                    (levenshtein cons-lines x) (levenshtein s1 s2) 
                    (catch-all (relative-entropy Prcons-lines Prx)) (catch-all (relative-entropy Prs1 Prs2))
                    (jensen-shannon Prcons-lines Prx) (jensen-shannon Prs1 Prs2)
                    ((shell/sh "perl" "/home/kitia/bin/gaisr/src/rnadist.pl") :out)])))))

;;looks over dataset with vector [(count cons) (count new cons) (count
;;s1) (count s2) (original lev) (s.waterman lev) (original KL)
;;(s.waterman KL) (original JSD) (s.waterman JSD) (bp-dist)] located in
;;temp.[neg|pos].clj which can be read using slurp. data represents
;;how different the individual structures from the alifold predicted
;;structure is to the RFAM consensus structure. slurped files made
;;from (io/with-out-writer "/home/kitia/bin/gaisr/src/temp.neg.clj"
;;(pr (vec (map vec negans))))
(map (fn [pos neg]
       (let [f (fn [i] (map (fn [[l _ _ len origlev s1s2lev _ _ origjsd s1s2jsd bpdist]]
                             [s1s2jsd (/ s1s2lev len)  origjsd (/ (Integer/parseInt bpdist) l)]) i))
             plev (f pos)
             nlev (f neg)]
         [(stats/mean (map first plev))
          (stats/mean (map first nlev))
          (stats/sd (map first plev))
          (stats/sd (map first nlev))
          (stats/mean (map second plev))
          (stats/mean (map second nlev))
          (stats/mean (map third plev))
          (stats/mean (map third nlev))
          (stats/mean (map #(->> % (drop 3) first) plev))
          (stats/mean (map #(->> % (drop 3) first) nlev))
          ]))
     (read-string (slurp "/home/kitia/bin/gaisr/src/temp.pos.clj")) (read-string (slurp "/home/kitia/bin/gaisr/src/temp.neg.clj")))

;;takes an input file and takes a percentage out for making a CM and
;;leaving the rest of testing.
(let [f "/home/kitia/bin/gaisr/alignments/ECL10_102611sort.sto"
              [_ seq-lines cons-lines] (join-sto-fasta-lines f "")
              n (count seq-lines)
              [a & b] (partition-all (max 3 (* n 0.3)) (first (repeatedly 1000 #(shuffle seq-lines))))]
              (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/temp.sto"
                                  (println "# STOCKHOLM 1.0\n")
                                  (doseq [[nm [_ sq]] a]
                                         (cl-format true "~A~40T~A~%" nm sq))
                                  (cl-format true "~A~40T~A~%" (ffirst cons-lines) (-> cons-lines first second second))
                                  (println "//"))
              (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/test.fa"
                (doseq [[nm [_ sq]] (apply concat b)]
                  (let [s (str/replace-re #"(\.|-)" "" sq)]
                  (println (str ">" nm))
                  (println s)
                  (println (str ">" nm "-bad"))
                  (println (markov-chain-seq s))))))












;;test sequences to test neutrality calc. stored in a file for reading
;;from later
(let [seqs (partition-all 2 (io/read-lines "/home/kitia/bin/gaisr/trainset2/test.fa"))]
  (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/neutrality.clj"
    (prn
     (vec (apply concat
                 (for [iseq (take 1 (partition-all 10 seqs))]
                   (vec
                    (pmap (fn [[name s] ]
                            (let [s (.toUpperCase s)
                                  nm (stats/mean (neutrality s))
                                  nc (neutrality_rand (fold s))]
                              [name nm (count (filter #(< nm %) nc)) nc]))
                          iseq))))))))

(defn markov-step
  "probs are a map with state and probability of occuring {:A 0.1 :H
  0.5 :T 0.2 :X 0.2}"
  
  [probs]
  (let [cum-probs (reductions + (vals probs))
        prob-table (map vector cum-probs (keys probs))
        r (rand)
        step (->> prob-table
                  (drop-while #(>= r (first %)))
                  first
                  second)]
    #_(prn prob-table) step))

(def mc
  {:states ["healthy" "sick"]
   :observation ["normal" "cold" "dizzy"]
   :start_probability {"healthy" 0.6 "sick" 0.4}
   :transition-probability {"healthy" {"healthy" 0.7 "sick" 0.3}
                            "sick" {"healthy" 0.4 "sick" 0.6}}
   :emission-probability {"healthy" {"normal" 0.5 "cold" 0.4 "dizzy" 0.1}
                          "sick" {"normal" 0.1 "cold" 0.3 "dizzy" 0.6}}})
(defn sim-mc [n]
  (let [start-state (let [start (->> (get mc :start_probability) markov-step)]
                      [start
                       (markov-step (get-in mc [:emission-probability start]))])
        next-state (fn [[cur-state emit]]
                     (let [trans (get-in mc [:transition-probability cur-state])
                           emis (get-in mc [:emission-probability cur-state])]
                       [(markov-step trans) (markov-step emis)]))
        ]
    (->> (iterate #(next-state %) start-state)
         (take n ))))
(defn viterbi
  "finds the most probable state path given a set of observations. for
  a test case let emissions = [\"normal\" \"cold\" \"dizzy\"]"

  [emissions]
  (let [hidden-states (keys (mc :start_probability))
                                        ;ans (sim-mc 3)
        y emissions ;["normal" "cold" "dizzy"] ;mc steps
        argmax (fn [x] (last (sort-by second x)))
        v1k (fn [k obs] (* (get-in mc [:emission-probability k obs]) ;pr(obs|init)
                          (get-in mc [:start_probability k]))) ;pr(init state)
        start (argmax (map (fn [state]
                             [state (v1k state (first y))])
                           hidden-states))]
    (loop [y (rest y)
           path (conj [] start)]
                                        ;(prn :path path :y y)
      (if (seq y)
        (let [obs (first y)
              vtk (fn [x k obs]
                    (* (get-in mc [:emission-probability k obs]);pr(obs|new)
                       (get-in mc [:transition-probability x k]);pr(old->new)
                       (-> path last second))) ;pr(old state)
              ]
          (recur (rest y)
                 (conj path
                       (argmax 
                        (map (fn [newstate]
                               [newstate (vtk (-> path last first) newstate obs)])
                             hidden-states)))))
        path)) ))
;;(pr old state * pr old->new * pr obs|new)

  

(def foo ;;takes ~2-2.5hrs to finish at 50 seqs 10 cpus
        (future (doall
           (pxmap
            (fn [i] (doall (map #(second (suboptimals % 10000)) i)))
            10
            (map mutant-neighbor
                 (map second (take 50 (take-nth 2 (partition-all 2 (io/read-lines "/home/peis/bin/gaisr/trainset2/test.fa"))))))))))

(let [wts (map second (take-nth 2 (partition-all 2 (io/read-lines "/home/peis/bin/gaisr/trainset2/test.fa"))))]
  (map (fn [wt mut]
         (map #(* (Math/sqrt (count wt))
                  (- 1 (pearsonsCC wt %))) mut))
       (map #(second (suboptimals % 10000)) wts) @foo))


;;generate 100 dinucleotide shuffles of the given sequence (from L10)
;;and then find the pSDC for each of these wt and shuffled versions
(def bar
  (future
    (doall
     (let [s (.toUpperCase "UCUAAAAGAACUGACCGAAGACAGUAGGGGACGAAAGUCAUAAACUUCCUACCgAGGACaAAUAUCAAAAUGAUA")
           [stc stv] (suboptimals s 10000)
           wts  (inverse-fold st 5)
           muts (pxmap
                 (fn [i] (map #(second (suboptimals % 10000)) i))
                 10
                 (map #(mutant-neighbor %) (cons s wts)))]
       (map (fn [wt mut]
              (map #(* (Math/sqrt (count wt))
                       (- 1 (pearsonsCC wt %))) mut))
            (cons stv (map #(second (suboptimals % 10000)) wts)) muts)))))

(def bar
  (future
    (doall
     (let [generate-vectors (fn [s]
                              (let [s (.toUpperCase s)
                                    [stc stv] (suboptimals s 10000)
                                    wts  (inverse-fold stc 100)
                                    muts (pxmap
                                          (fn [i] (doall (map #(second (suboptimals % 10000)) i)))
                                          10
                                          (map #(mutant-neighbor %) (cons s wts)))] 
                                [(cons stv (map #(second (suboptimals % 10000)) wts)) muts]))
           inseqs (let [fdir "/home/peis/bin/gaisr/trainset2/"
                        fvector ["RF00555-seed.4.fasta" "RF00558-seed.3.fasta"
                                 "RF00559-seed.7.fasta" "RF00167-seed.9.fasta"]
                        randseq (fn [coll] (first (last (repeatedly 100 #(shuffle coll)))))]
                    (conj (map (fn [f]
                                 (->> (io/read-lines (str fdir f))
                                      rest
                                      (take-nth 2)
                                      ))
                               fvector) 
                          '("UCCGGAAUAUCCUGUCCGAGAUUGUGGGUGAUACuGUUUgaGUAUCUUAAuCAAAAAaaCCUGCAUgAGACUGGGGUAAGAACUGUAG" "UCUAAAAGAACUGACCGAAGACAGUAGGGGACGAAAGUCAUAAACUUCCUACCgAGGACaAAUAUCAAAAUGAUA")))
           ]
       (map (fn [k vs]
              (assoc {} k (for [i (map generate-vectors vs)]
                            (map psdc i))))
            [:L10 :L13 :L20 :L21 :FMN] inseqs)))))





;;load up database stuff for usage (use '[edu.bc.bio.gaisr.db-actions :only [sql-query mysql-ds]])
;; (sql-query
;;            "select be.name,an.ancestors
;;           from bioentry as be, taxon as tx, ancestor as an
;;           where be.taxon_id=tx.taxon_id
;;           and   tx.ncbi_taxon_id=an.ncbi_taxon_id
;;           and   be.name = \"NC_009654\"")


;;make graphs 
(let [abc (map #(map (fn [x] (stats/mean x)) (partition-all 3 %)) barr)
      wt (first abc)
      n (count wt)
      xy (charts/xy-plot (range 1 (inc n)) wt 
                         :y-label "pSDC"
                         :x-label "position"
                         :title "seq3 collapsed mutations" 
                         :legend true 
                         :series-label "wt"
                         :points true)]
  (view xy) 
  (map (fn [i y]
         (charts/add-lines xy (range 1 (inc n)) y :series-label (str "neg" i)))
       (iterate inc 1) (rest abc)))


;;create charts using highcharts

(defn most-likely-ancestor
  "take a vector of strings and finds the consensus sequence by
  choosing the most common base at a position."

(let [s "AUGCUAGUACGGAU"
      n (count s)
      ztable (atom {})
      Z (fn Z [i j S]
          (let [S (.toUpperCase S)
                n (- j i -1)
                bp? (fn [b1 b2] ;;b1=base1 b2=base2
                      (let [bp #{"AU" "UA" "GC"  "CG" "GU" "UG"}]
                        (contains? bp (str b1 b2))))
                e (fn [i j S] (let [s1 (subs S i (inc i))
                                   s2 (subs S j (inc j))
                                   R 2
                                   T 310
                                   E (fn [b1 b2] (if (bp? b1 b2) 1 0))] ;;Energy of basepair, E(basepair)
                               (E s1 s2) #_(Math/exp (/ (E s1 s2) R T -1))))
                u 0 ;;min loop size
                result (if (<= (- j i) u) 1
                           (+ (get @ztable [i (dec j)] (Z i (dec j) S)) ;;j unpaired
                              (* (e i j S) (get @ztable [(inc i) (dec j)] (Z (inc i) (dec j) S))) ;;i,j pair
                              (reduce (fn [x k]  ;;k,j paired for an intermediate a<k<b
                                        (+ x (* (e k j S) (get @ztable [i (dec k)] (Z i (dec k) S)) (get @ztable [(inc k) (dec j)] (Z (inc k) (dec j) S)))))
                                      0 (range (inc i) (- j u)))))]
            (swap! ztable #(assoc % [i j] result)) result))]
  #_(for [i (range 1 n)
          j (range n)
          :when (< (+ j i) n)]
    (swap! ztable #(assoc % [j (+ i j)] (Z j (+ i j) s))))
  (Z 0 5 s) @ztable)




(deref (future #((fn [n target] (sh "RNAinverse"
                                   "-Fmp"
                                   (str "-R" n)
                                   "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                   :in target)) 2 "((...))")) 10 :error)


(defn future-timeout [f timeout-in-s]
  (.get f timeout-in-s (java.util.concurrent.TimeUnit/SECONDS)))

(.get (future #((fn [n target] (sh "RNAinverse"
                                            "-Fmp"
                                            (str "-R" n)
                                            "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                            :in target)) 2 "((...))")) 10 (java.util.concurrent.TimeUnit/SECONDS))

(deref (future-call #((fn [n target] (sh "RNAinverse"
                                        "-Fmp"
                                        (str "-R" n)
                                        "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par"
                                        :in target)) 2 "((...))")) 1000 :error)

(defmacro find-string [s1 s2]
  '(re-find s1 ~s2))

(let [known-primes #{}
      some-set (range 2 100000)]
  (loop [i some-set
         kp known-primes]
    (letfn [(prime? [n]
              (or (contains? kp n)
                  (every? #(pos? (mod n %)) (range 2 (int (Math/sqrt n))))))]
      (if (seq i)
        (recur (rest i)
               (if (prime? (first i)) (conj kp (first i)) kp))
        kp))))
