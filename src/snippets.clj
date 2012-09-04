(ns snippet
  (:require [clojure.contrib.string :as str]
            [clojure.java.shell :as shell]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs]
            [incanter.stats :as stats]
            [clojure.set :as sets])
  (:use edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.utils
        edu.bc.utils.probs-stats
        ;smith_waterman
        [clojure.contrib.pprint
         :only (cl-format compile-format)]))

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
   random sequence th at shares a thr of sequence identity"

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

(defn svm-features [f]
  (let [m (profile (read-sto f))
        mi (fn [x]
             ((group-by (fn [[[i j] _]]
                          (contains? (set
                                      (apply concat (map #(vec %) (m :pairs))))
                                     [i j]))
                        x) true))]
    (prn "zscore" (stats/mean (zscore m)))
    (prn "sci" (sci (energy-of-aliseq2 m) (energy-of-seq2 m)))
    (prn "information" (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))))
    (prn "mutual info" (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)))
    (prn "pairwise identity" (let [id (pairwise_identity (m :seqs))]
                               (double (/ (apply + (map first id))
                                          (apply + (map second id))))))
    (prn "number of seqs" (count (m :seqs)))))


(do (io/with-out-writer "/home/kitia/bin/gaisr/trainset2/train2.csv"
      (println "zscore, sci, information, MI, JS, pairwise identity, number of seqs, class")
      (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/pos/list.txt")] 
        (let [m (profile (read-sto (str "/home/kitia/bin/gaisr/trainset2/pos/" f)))
              mi (fn [x]
                   ((group-by (fn [[[i j] _]]
                                (contains? (set
                                            (apply concat (map #(vec %) (m :pairs))))
                                           [i j]))
                              x) true))]
          (print (stats/mean (zscore m)) ",")
          (print (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ",")
          (print (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))) ",")
          (print (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)) ",")
          (print (stats/mean (information_only_bp (reduce (fn [x i]
                                                          (assoc x i (gutell_calcs/JS (col->prob (->> (transpose (m :seqs))
                                                                                         (drop i)
                                                                                         first)
                                                                                    :gaps true) :Q (m :background))))
                                                        {} (range (m :length))) (m :pairs))) ",")
          (print (let [id (pairwise_identity (m :seqs))]
                   (double (/ (apply + (map first id))
                              (apply + (map second id))))) ",")
          (print (count (m :seqs)) ",1\n")))
    (doseq [f (io/read-lines "/home/kitia/bin/gaisr/trainset2/neg/list.txt")] 
      (let [dir "/home/kitia/bin/gaisr/trainset2/neg/"
            sto (str (subs f 0 (- (count f) 3)) "sto")
            m (profile (read-sto (snippet/aln->sto (str dir f) (str dir sto))))
            mi (fn [x]
                 ((group-by (fn [[[i j] _]]
                              (contains? (set
                                          (apply concat (map #(vec %) (m :pairs))))
                                         [i j]))
                            x) true))]
        (print (stats/mean (zscore m)) ",")
        (print (sci (energy-of-aliseq2 m) (energy-of-seq2 m)) ",")
        (print (stats/mean (information_only_bp (gutell_calcs/entropy m) (m :pairs))) ",")
        (print (stats/mean (->> (mi (gutell_calcs/mutual_info m)) (into {}) vals)) ",")
        (print (stats/mean (information_only_bp (reduce (fn [x i]
                                                          (assoc x i (gutell_calcs/JS (col->prob (->> (transpose (m :seqs))
                                                                                         (drop i)
                                                                                         first)
                                                                                    :gaps true) :Q (m :background))))
                                                        {} (range (m :length))) (m :pairs))) ",")
        (print (let [id (pairwise_identity (m :seqs))]
                 (double (/ (apply + (map first id))
                            (apply + (map second id))))) ",")
        (print (count (m :seqs)) ",0\n")))))

(defn aln->sto
  "takes an alignment in Clustal W format and produces a sto file by using RNAalifold to determine
   the structure and then making it into a sto file adding header and a consensus line"

  [aln sto & {fold_alg :fold_alg :or {Q "RNAalifold" }}]
  (if (= fold_alg "RNAalifold")
    (let [st (->> ((shell/sh "RNAalifold"  "-P" "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par" "-r" "--noPS" aln) :out)
                  (str/split-lines)
                  second
                  (str/split #" ")
                  first)
          sq (rest (second (join-sto-fasta-lines aln "")))]
      (io/with-out-writer sto
        (println "# STOCKHOLM 1.0\n")
        (doseq [[n [_ s]] sq]
          (cl-format true "~A~40T~A~%" n (str/replace-re #"\-" "." s)))
        (cl-format true "~A~40T~A~%" "#=GC SS_cons" st)
        (println "//"))
      sto)
    (shell/sh "perl" "/home/kitia/bin/gaisr/src/mod_cmfinder.pl" aln sto)))

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

(defn randseqs
  "takes seq-lines and the number of random sequences to draw from it"
  [seq-lines n]
  (let [rand-sqs (map (fn [x]
                        (nth seq-lines x))
                      (for [i (range n)]
                        (rand-int (count seq-lines))))]
    (map (fn [[nm [_ sq]]]
           [nm sq])
         (sort-by #(-> % second first) <
                  (map (fn [[nm [uid sq]] ungap-sq]
                         [nm [uid ungap-sq]])
                       rand-sqs (->> (map (fn [[nm [uid sq]]]
                                            sq)
                                          rand-sqs)
                                     transpose
                                     (remove #(empty? (str/replace-re #"\.|\-" "" %)))
                                     transpose))))))

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

(defn sto->fasta
  "converts a sto file to a fasta file. removes the gaps from the
   sequences. This should be used in order to make a fasta file for
   input into CMfinder"
  
  ([sto & {type :type :or {type :sto }}]
     (let [lines (second (join-sto-fasta-lines sto ""))]
       (doseq [[nm [_ sq]] lines]
         (println (str ">" nm))
         (println (str/replace-re #"\." "" sq)))))
  
  ([sto outfasta]
     (io/with-out-writer outfasta
      (sto->fasta sto))))

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

(defn fold
  "Folds a sequence of RNA and returns only the target structure"
  
  [s]
  (->> ((shell/sh "RNAfold"  "-P" "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par" "--noPS" :in s) :out)
                         (str/split-lines)
                         second
                         (str/split #" ")
                         first))
(defn neutrality
  "takes a string s and finds the 3L mutants. String can only contain
   letters A, C, G, U. returns neutrality of each of the seq compared
   to each of the 3L 1 neighbors mutants"

  [s]
  (let [b {"A" 1 "G" 1 "U" 1 "C" 1}
        st (fold s)
        neighbors (fn [s]
                    (flatten
                     (for [i (range (count s))]
                       (map (fn [r]
                              (str (subs s 0 i) r (subs s (inc i) (count s))))
                            (keys (dissoc b (subs s i (inc i)))))))) ]
    (map (fn [neighbor]
           (/ (- (count s) (levenshtein st (fold neighbor)))
              (count s)))
         (neighbors s))))

(defn neutrality_rand
  "neutrality from finding a random seq with the same structure. if neutrality
   for native seq > rand seq then the native seq is more robust."
  
  ([target]
     (neutrality_rand target 100))
  
  ([target n]
     (let [cand (loop [c 0
                       cand []]
                  (if (< c n)
                    (let [x (remove nil?
                                    (flatten 
                                     (map (fn [[s ensemble]]
                                            (when-not (re-find #"d=" s) (re-find #"\w+" s)))
                                          (partition-all 2 (str/split-lines ((shell/sh "RNAinverse" "-Fmp" (str "-R" (- n c)) "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par" :in target) :out))))))]
                      (recur (count (distinct cand)) (concat (distinct cand) x)))
                    (distinct cand)))]
      
       #_(prn distinct? cand)
       ;;(filter #(= target (fold %)) cand)
       (map #(stats/mean (neutrality %)) cand))))


;;est sequences to test neutrality calc. stored in a file for reading
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


