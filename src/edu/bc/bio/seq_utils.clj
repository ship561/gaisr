;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             S E Q - U T I L S                            ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011-2012 Trustees of Boston College                       ;;
;;                                                                          ;;
;; Permission is hereby granted, free of charge, to any person obtaining    ;;
;; a copy of this software and associated documentation files (the          ;;
;; "Software"), to deal in the Software without restriction, including      ;;
;; without limitation the rights to use, copy, modify, merge, publish,      ;;
;; distribute, sublicense, and/or sell copies of the Software, and to       ;;
;; permit persons to whom the Software is furnished to do so, subject to    ;;
;; the following conditions:                                                ;;
;;                                                                          ;;
;; The above copyright notice and this permission notice shall be           ;;
;; included in all copies or substantial portions of the Software.          ;;
;;                                                                          ;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          ;;
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       ;;
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    ;;
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   ;;
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ;;
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    ;;
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          ;;
;;                                                                          ;;
;; Author: Jon Anthony                                                      ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.bio.seq-utils

  "Intended to include/contain various ancillary operations on
   sequence data. Refactoring will sooner than later make this
   obsolete.  Note it now contains what once was in seq-utils2 and
   that should be factored as well."

  (:require clojure.string
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.contrib.math :as math]
            [edu.bc.fs :as fs])

  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.stats :only (hamming-distance)]
        ))






(defn sto-GC-and-seq-lines [stofilespec]
  (seq/separate
   #(and (> (count %) 1)
         (or (not (.startsWith % "#"))
             (or (.startsWith % "#=GC SS_cons")
                 (.startsWith % "#=GC RF"))))
   (io/read-lines (fs/fullpath stofilespec))))


(defn join-sto-fasta-lines [infilespec origin]
  (let [[seqcons-lines gc-lines] (sto-GC-and-seq-lines infilespec)
        gc-lines (if (not= origin "")
                   (concat (take 1 gc-lines) [origin] (drop 1 gc-lines))
                   gc-lines)
        recombined-lines (sort-by
                          #(-> % second first)
                          (vec (reduce
                                (fn [m l]
                                  (let [[nm sq]
                                        (cond
                                         ;;splits the line apart and
                                         ;;hopefully creates vector
                                         ;;["#GC SS_cons" structure]
                                         (.startsWith l "#=GC SS_cons")
                                         [(str/join " " (butlast (str/split #"\s+" l)))
                                          (last (str/split #"\s+" l))]

                                         (.startsWith l "#")
                                         (str/split #"\s{2,}+" l)

                                         :else
                                         (str/split #"\s+" l))

                                        prev (get m nm [(gen-uid) ""])]
                                    (assoc m  nm [(first prev)
                                                  (str (second prev) sq)])))
                                {} seqcons-lines)))
        {seq-lines false cons-lines true} (group-by
                                           #(or (.startsWith (first %) "//")
                                                (.startsWith (first %) "#"))
                                           recombined-lines)]
    [gc-lines seq-lines cons-lines]))


(defn join-sto-fasta-file
  "Block/join unblocked sequence lines in a sto or fasta file. For sto
   files ORIGIN is a #=GF line indicating tool origin of file.  For
   example, '#=GF AU Infernal 1.0.2'. Defaults to nothing."

  [in-filespec out-filespec
   & {origin :origin :or {origin ""}}]
  (let [[gc-lines seq-lines cons-lines]
        (join-sto-fasta-lines in-filespec origin)]
    (io/with-out-writer (fs/fullpath out-filespec)
      (doseq [gcl gc-lines] (println gcl))
      (doseq [sl seq-lines]
        (let [[nm [_ sq]] sl]
          (cl-format true "~A~40T~A~%" nm sq)))
      (doseq [cl cons-lines]
        (let [[nm [_ sq]] cl]
          (cl-format true "~A~40T~A~%" nm sq))))))

(defn nsubseq-sto-file
  "Reduces the number of sequences in an unblocked sto file. Takes in
   a sto file and outputs a new sto file with nseq number of
   lines. The function currently uses hamming distance to determine
   the difference between any 2 sequnces in the alignment. After the
   pairwise hamming distances are determined, the alignment is then
   re-ordered according to these hamming distance vectors. These
   vectors should allow similar sequences to be grouped together. Once
   grouped together, the function goes through and eliminates a
   sequence when comparing seq i and i+1. i+1 is removed if the two
   sequences are less different than a set threshold, thr=0.5
   initially. The loop compares i and i+1 then i+2 vs i+3 and so
   forth. The threshold will increase by 0.01 each loop so that the
   final list of sequences will be nseq long."

  [infile outfile nseq]
  (let [[gc-lines seq-lines cons-lines] (join-sto-fasta-lines infile "")
        ;;create a map where k=name v=[hamming distance vector seq]
        calc-diff (fn [inseqs] (reduce (fn [m [n [_ s1]]]
                                   (assoc m n
                                          [(vec (for [[_ [_ s2]] inseqs]
                                                  (hamming-distance s1 s2)))
                                           s1]))
                                 {} inseqs))
        ;;returns a map k=name v=seq. removes the hamming distance
        ;;from the vector
        remove-diff (fn [x] (into {} (map (fn [[nm [_ sq]]]
                                           [nm sq])
                                         x)))
        ;;output of new sto file to outfile
        print-sto (fn [x]
                    (io/with-out-writer outfile
                      (doseq [gc gc-lines] (println gc))
                      (doseq [sl x]
                        (let [[nm [_ sq]] sl]
                          (cl-format true "~A~40T~A~%" nm sq)))
                      (doseq [cl cons-lines]
                        (let [[nm [_ sq]] cl]
                          (cl-format true "~A~40T~A~%" nm sq)))))]
    (print-sto
      (loop [m (into {} (sort-by #(first (last %)) (calc-diff seq-lines)))
             nm m
             thr 0.05]
        (let [i (first (keys m))
              j (second (keys m))
              [_ s1] (m i)
              [_ s2] (m j)
              diff (/ (hamming-distance s1 s2) (count s1))]
          (if (> (count nm) nseq)
            (if (> (count m) 3)
              (recur (dissoc m i j)
                     (if (< diff thr) (dissoc nm j) nm)
                     thr)
              (recur nm
                     (if (< diff thr) (dissoc nm j) nm)
                     (+ thr 0.01)))
            nm))))))

 
(defn check-sto
  "Checks a sto file to ensure that there are valid characters being
   used in the sequences consensus structure line. Will print out
   errors in the sto file by sequence number.  Input requires a sto
   file"
  
  [sto]
  (let [valid-symbols #{"A" "C" "G" "U"
                        "-" "." ":" "_"
                        "a" "b" "B"
                        "(" ")" "<" ">"}
        [_ seq-lines cons-lines] (join-sto-fasta-lines sto "")
        [_ [_ cl]] (first cons-lines)
        ;;adds numbers to the sequences so that they
        ;;can be identified
        sl (partition 2 
                      (interleave (iterate inc 1)
                                  (reduce (fn [v [_ [_ sq]]]
                                            (conj v (.toUpperCase sq)))
                                          [] seq-lines)))
        ;;checks for valid symbols
        check-char (fn [[n s]] 
                     [n (every? #(contains? valid-symbols %) (rest (str/split #"" s)))])
        ;;checks to see if sequences have same
        ;;length as consensus structure
        check-len (fn [[n s]] 
                    [n (= (count s) (count cl))])]
    (cond
     (some #(false? (second %)) (map #(check-char %) sl))
     (prn "sequence contains invalid character in: " (remove #(second %) (map #(check-char %) sl)))
     
     (some #(false? (second %)) (map #(check-len %) sl))
     (prn "sequence contains invalid length compared to cons-line in: " (remove #(second %) (map #(check-len %) sl)))
     
     (false? (second (check-char [1 cl])))
     (prn "consensus structure contains invalid character")
     
     :else
     (prn "sto is good"))))
        

;;; Convert STO format to ALN format (ClustalW format).  This is
;;; needed by some processors which cannot take a Stockholm alignment
;;; format but need an "equivalent" in ClustalW ALigNment format.
;;;
;;; OK, (9-Feb-2012) some tools seem to need things blocked while
;;; others don't work if they are blocked.  Worse, what counts as
;;; valid Clustal/aln format or not is ill defined with multiple
;;; definitions in the community (e.g., many claim 60 character seqs
;;; per line but others say 50; some claim must be blocked, others say
;;; unblocked is valid).  So, we have two variants.  One which blocks
;;; and a main driver which calls blocked version if requested or just
;;; does simple unblocked itself.
;;;
(defn sto->aln-blocked
  "Convert a stockhom format alignment file into its ClustalW
   equivalent BLOCKED ALN format. Blocking is done in 60 character
   chunks.  STOIN is the filespec for the stockholm format file and
   ALNOUT is the filespec for the resulting conversion (it is
   overwritten if it already exists!)"

  [stoin alnout]
  (let [seq-lines (second (join-sto-fasta-lines stoin ""))
        seq-lines (map (fn [[nm [uid sl]]]
                         [nm [uid (partition-stg
                                   60 (str/replace-re #"\." "-" sl))]])
                       seq-lines)]
    (io/with-out-writer alnout
      (println "CLUSTAL W (1.83) multiple sequence alignment\n")
      (loop [x seq-lines]
        (let [[nm [uid sl]] (first x)]
          (when (not-empty sl)
            (do
              (doseq [[nm [uid sl]] x]
                  (cl-format true "~A~40T~A~%" nm (first sl)))
              (println "")
              (recur (map (fn [[nm [uid sl]]]
                            [nm [uid (rest sl)]])
                          x)))))))
    alnout))


(defn sto->aln
  "Convert a stockhom format alignment file into its ClustalW
   equivalent ALN format.  STOIN is the filespec for the stockholm
   format file and ALNOUT is the filespec for the resulting
   conversion (it is overwritten if it already exists!)

   BLOCKED is a boolean indicating whether the output should be
   blocked (60 chars per chunk).  Default is unblocked."

  [stoin alnout & {blocked :blocked :or {blocked false}}]
  (if blocked
    (sto->aln-blocked stoin alnout)
    (let [seq-lines (filter #(not (or (.startsWith % "//") (re-find #"^#" %)))
                            (first (sto-GC-and-seq-lines stoin)))
          seq-lines (map #(str/replace-re #"\." "-" %) seq-lines)]
      (io/with-out-writer (fs/fullpath alnout)
        (println "CLUSTAL W (1.83) multiple sequence alignment\n")
        (doseq [sl seq-lines]
          (println sl)))
      alnout)))

(defn reverse-compliment
  "Reverse compliment of DNA/RNA sequence SQ. Return the sequence that
   would pair with (reverse SQ).

;;; ----------------------------------------------------------------------
;;;
;;; Convert Sto and Fasta split sequence format files into conjoined
;;; versions.  Many Sto and Fasta files from various sites come in old
;;; fashioned 80 col mode where sequences are split at 80 column mark.
;;; For aligned files this is even worse as you have groups of
;;; sequences split across lines separated by whole pages of other
;;; sequence (parts).  For example, RFAM alignments.  This group puts
;;; all those back together so that each sequence is on a single line.

   Ex: (reverse-compliment \"AAGGAAUUCC\") => \"GGAAUUCCUU\"
  "
  [sq]
  (let [ucp (first (str/codepoints "U"))
        m (if (neg? (.indexOf sq ucp))
            {\A \T \T \A \G \C \C \G}
            {\A \U \U \A \G \C \C \G})]
    (str/join "" (reverse (map m sq)))))


;;; ----------------------------------------------------------------------
;;;
;;; Sequence shuffling with nucleotide/aa frequency distributions
;;; preserved.  Currently only provides support for dinucleotide
;;; distributions.  Which is basically bogus, but the most typical use
;;; case??


(defn nt-cntn
  "Frequencies of N-(nucleotides/amino-acids) of sequence SEQUENCE,
   i.e., for N=2 and bases freq of dinucleotides.  Effectively a
   rename of utils:freqn.
   "
  [n sequence]
  (freqn n sequence))


(def test-nts
     (nt-cntn 2 (str/join "" (repeat 1 "acagtcaacctggagcctggt"))))

(def nts-map
     (reduce #(assoc %1 (subs (first %2) 0 1)
                     (conj (get %1 (subs (first %2) 0 1) []) %2))
             {} test-nts))

(def legal-ends
     (keep #(when (= (subs (first %1) 1 2) "t") %1) test-nts))


(def limit-counter (atom 50000))
(def print-limit-info (atom nil))

(defn check-limit [start len ss newseq nts-map]
  (swap! limit-counter dec)
  (when (= limit-counter 0)
    (when print-limit-info
      (prn start len ss newseq nts-map))
    (raise :type :explode)))


(def shuffles nil)

(defn dint-shuffle [start len nts-map legal-ends ss newseq & base-starters]
  (cond
   (and (= len 0) (some #(= (last newseq) (first %1)) legal-ends))
   (str ss (subs (last newseq) 1 2))

   (= len 0) :fail

   (nil? (nts-map start)) :fail

   :otherwise
   (loop [starters (or (first base-starters) (nts-map start))
          [pick cnt] (first starters)]
     (let [newstarters (drop 1 starters)
           mod-edges (keep  #(if (not= pick (first %1)) %1
                               (let [k (first %1) cnt (dec (second %1))]
                                 (when (> cnt 0) [k cnt])))
                            (nts-map start))
           new-map (if (not-empty mod-edges)
                     (assoc nts-map start
                            (concat (drop 1 mod-edges)
                                    (take 1 mod-edges)))
                     (dissoc nts-map start))
           result (dint-shuffle (subs pick 1 2) (dec len) new-map
                                legal-ends (str ss start)
                                (conj newseq pick))]
       (if (and result (not= result :fail)
                (not (@shuffles result)))
         result
         (do
           (check-limit start len ss newseq nts-map)
           (if (not-empty newstarters)
             (recur newstarters
                    (first newstarters))
             :fail)))))))


(defn dint-seq-shuffle [limit s]
  (let [len (count s)
        start (subs s 0 1)
        end (subs s (dec len) len)
        di-nts (nt-cntn 2 s)
        nts-map (reduce #(assoc %1 (subs (first %2) 0 1)
                                (conj (get %1 (subs (first %2) 0 1) []) %2))
                        {} di-nts)
        legal-ends (keep #(when (= (subs (first %1) 1 2) end) %1) di-nts)]
    (binding [shuffles (atom {})
              limit-counter limit-counter]
      (handler-case :type
        (loop [_ limit
               base-starters (nts-map start)]
          (let [s (dint-shuffle
                   start (dec len) nts-map legal-ends "" [] base-starters)]
            (when (not= :fail s)
              (swap! shuffles
                     assoc s (inc (get shuffles s 0))))
            (when (> _ 1)
              (recur (dec _)
                     (concat (drop 1 base-starters)
                             (take 1 base-starters))))))
        (handle :explode
          (println :limit-reached)))
        @shuffles)))




;;; Simple Monte Carlo subseq frequency approximation within total seqs
;;;

(defn alphabet [type]
  {:pre [(in type [:dna :rna :amino])]}
  (case type
        :dna ["A" "T" "G" "C"]
        :rna ["A" "U" "G" "C"]
        :amino ["A" "R" "N" "D" "C" "Q" "E" "G" "H" "I"
                "L" "K" "M" "F" "P" "S" "T" "W" "Y" "V"]))


(defn gen-random-bioseq  [type size & seed]
  {:pre [(> size 0)]}
  (let [alphabet (alphabet type)
        seed (first seed)
        seed (if seed (vec (str/upper-case seed)) [])
        size (- size (count seed))]
    (loop [s seed
           cnt 0]
      (if (< cnt size)
        (recur (conj s (rand-nth alphabet))
               (inc cnt))
        [(str/join "" s) s]))))


(defn monte-carlo-subseq-cnts-body [type sample-size seqsize sseq len seed]
  (loop [ss sample-size
         cnts {}
         sseqcnts {}]
    (if (= ss 0)
      [cnts sseqcnts]
      (let [s (first (gen-random-bioseq type seqsize (when seed sseq)))
            scnts (nt-cntn len s)]
        (recur (dec ss)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       cnts scnts)
               (assoc sseqcnts (scnts sseq) (inc (get scnts sseq 0))))))))


(defnk monte-carlo-subseq-cnts
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil]
  {:pre [(> size 0)
         (let [abet (alphabet type)]
           (every? #(in (str %1) abet) (str/upper-case sseq)))]}
  (let [len (count sseq)
        sseq (str/upper-case sseq)]
    (monte-carlo-subseq-cnts-body type sample-size size sseq len seed)))


(defn reduce-seqcnt-results [cntmaps]
  (loop [cntmaps cntmaps
         total-cnts {}
         sseqcnts {}]
    (if (empty? cntmaps)
      [total-cnts sseqcnts]
      (let [[tcnts scnts] (first cntmaps)]
        (recur (drop 1 cntmaps)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       total-cnts tcnts)
               (reduce (fn[m [k v]]
                         (assoc m k (+ (get m k 0) v)))
                       sseqcnts scnts))))))

(defnk pmonte-carlo-subseq-cnts
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil :par 10]
  {:pre [(> size 0)
         (and (integer? par) (<= 1 par 12))
         (let [abet (alphabet type)]
           (every? #(in (str %1) abet) (str/upper-case sseq)))]}
  (let [len (count sseq)
        sseq (str/upper-case sseq)
        psz (floor (/ sample-size par))
        [q r] (div sample-size psz)
        ssizes (repeat q psz)
        ssizes (conj (drop 1 ssizes) (+ (first ssizes) r))]
    (reduce-seqcnt-results
     (pmap #(monte-carlo-subseq-cnts-body type %1 size sseq len seed)
           ssizes))))


(defnk monte-carlo-subseq-freq
  [:type :dna :size 1024 :sseq "" :sample-size 10000 :seed nil :par 10]
  (let [f (if (and par (integer? par) (> par 1))
            pmonte-carlo-subseq-cnts
            monte-carlo-subseq-cnts)
        [cnts sseqcnts] (f :type type :size size :sseq sseq
                           :sample-size sample-size :seed seed :par par)
        sseq-total (cnts (str/upper-case sseq))
        total (reduce (fn[cnt [k v]] (+ cnt v)) 0 cnts)]
    [sseq-total total (double (/ sseq-total total)) sseqcnts]))



