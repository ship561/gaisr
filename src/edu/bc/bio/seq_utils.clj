;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                             S E Q - U T I L S                            ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011 Trustees of Boston College                            ;;
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

  "Intended to include/contain various utilities for operating with and on
   sequence data.  Originally this was all about working with (NCBI) Blast+,
   CMFINDER (and associated utils) and Infernal tools"

  (:require clojure.string
            [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clojure.xml :as xml]
            [clj-shell.shell :as sh]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        [incanter.stats :only (hamming-distance)]
        ))




;;; ----------------------------------------------------------------------
;;;
;;; Convert Sto and Fasta split sequence format files into conjoined
;;; versions.  Many Sto and Fasta files from various sites come in old
;;; fashioned 80 col mode where sequences are split at 80 column mark.
;;; For aligned files this is even worse as you have groups of
;;; sequences split across lines separated by whole pages of other
;;; sequence (parts).  For example, RFAM alignments.  This group puts
;;; all those back together so that each sequence is on a single line.


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
                     [n (every? #(contains? valid-symbols %)
                                (rest (str/split #"" s)))])
        ;;checks to see if sequences have same
        ;;length as consensus structure
        check-len (fn [[n s]]
                    [n (= (count s) (count cl))])]
    (cond
     (some #(false? (second %)) (map #(check-char %) sl))
     (prn "sequence contains invalid character in: "
          (remove #(second %) (map #(check-char %) sl)))

     (some #(false? (second %)) (map #(check-len %) sl))
     (prn "sequence contains invalid length compared to cons-line in: "
          (remove #(second %) (map #(check-len %) sl)))

     (false? (second (check-char [1 cl])))
     (prn "consensus structure contains invalid character")

     :else
     (prn "sto is good"))))




;;; ----------------------------------------------------------------------
;;; BLAST functions.  These are currently based on NCBI BLAST+ (which
;;; purportedly sort of sucks, but is still open source where wu-blast
;;; is now commercial via AB-Blast).

(def default-binary-db
     (or (getenv "MLAB_DEFAULT_BINARY_DB")
         "/data2/BioData/BlastDBs/other_genomic"))


(defn gen-entry-file [entries file]
  (io/with-out-writer (io/file-str file)
    (doseq [e entries]
      (println e)))
  file)


(defn has-loc? [entries]
  (let [entries (cond
                 (seq? entries) entries
                 (fs/exists? entries) (str/split #"\n" (slurp entries))
                 :else (raise :type :unknown-entries :entries entries))]
    (let [e (first (ensure-vec entries))]
      (re-find #"[0-9]+\-[0-9]+" e))))




(defn nms-sqs->fasta-file [nms-sqs filespec]
  (let [filespec (fs/fullpath filespec)]
    (io/with-out-writer filespec
      (doseq [[nm sq] nms-sqs]
        (println (str ">" nm))
        (println sq)))
    filespec))


(defn entry-range->fasta-file [entry range blastdb & [filespec]]
  (let [fasta-filespec (if filespec filespec (fs/tempfile entry ".fna"))
        blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")]
    (assert-tools-exist [blastdbcmd])
    (let [[range strand] (str/split #"/" range)
          strand (if (= "-1" strand) "minus" "plus") ; => default plus
          cmdargs ["-db" blastdb
                   "-entry" entry "-range" range "-strand" strand
                   "-line_length" "14000000"
                   "-outfmt" "%f" "-out" (fs/fullpath fasta-filespec)]]
      (runx blastdbcmd cmdargs)
      fasta-filespec)))


(defn- entry-file->fasta-file-ranges [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        tmp-file (fs/tempfile "fasta-out" ".fna")]
    (assert-tools-exist [blastdbcmd])
    (io/with-out-writer fasta-filespec
      (do-text-file [efile]
        (let [[entry range] (str/split #" " $line)
              [range strand] (str/split #"/" range)
              strand (if (= "-1" strand) "minus" "plus") ; => default plus
              cmdargs ["-db" blastdb
                       "-entry" entry "-range" range "-strand" strand
                       "-line_length" "14000000"
                       "-outfmt" "%f" "-out" tmp-file]]
          (catch-all (runx blastdbcmd cmdargs))
          (do-text-file [tmp-file]
                        (println $line)))))
    (fs/rm tmp-file)
    fasta-filespec))


(defn- entry-file->fasta-file-full [efile fasta-filespec blastdb]
  (let [blast-path (get-tool-path :ncbi)
        blastdbcmd (str blast-path "blastdbcmd")
        cmdargs ["-db" blastdb "-entry_batch" efile
                 "-line_length" "14000000"
                 "-outfmt" "%f" "-out" fasta-filespec]]
    (assert-tools-exist [blastdbcmd])
    (catch-all (runx blastdbcmd cmdargs))
    fasta-filespec))


(defnk entry-file->fasta-file [efile :loc nil :blastdb nil]
  (let [efile (fs/fullpath efile)
        filespec (fs/fullpath (fs/replace-type efile ".fna"))
        blastdb (if blastdb blastdb default-binary-db)]
    (if loc
      (entry-file->fasta-file-ranges efile filespec blastdb)
      (entry-file->fasta-file-full efile filespec blastdb))
    filespec))


(defnk entry-file->blastdb [efile :out nil :blastdb nil :loc nil :ids nil]
  (let [entries-fasta (entry-file->fasta-file efile :loc loc :blastdb blastdb)
        blast-db-file (fs/fullpath
                       (if out out (fs/replace-type efile ".blastdb")))
        blast-path (get-tool-path :ncbi)
        mkdbcmd (str blast-path "makeblastdb")
        cmdargs ["-parse_seqids"
                 "-in" entries-fasta "-dbtype" "nucl"
                 "-out" blast-db-file]
        cmdargs (if ids cmdargs (vec (rest cmdargs)))]
    (assert-tools-exist [mkdbcmd])
    (runx mkdbcmd cmdargs)
    blast-db-file))




(defn get-selection-fna
  "Get a fasta file for SELECTIONS.  If selections is a file, assumes
   it is in fact the fasta file.  If selections is a collection,
   assumes the collection is a set of pairs [nm sq], and converts to
   corresponding fasta file.  Returns full path of result file."
  [selections]
  (if (coll? selections)
    (nms-sqs->fasta-file selections)
    (fs/fullpath selections)))


(defnk blast
  "Implementation of blast requests.  See blastn, tblastn, etc. for details
   on main arguments.  PGM is a keyword denoting which blast program."
  [pgm in
   :out nil
   :blastdb nil
   :strand :plus
   :word-size 8
   :evalue 10
   :fmt "10 qseqid qstart qend evalue sseqid sstart send"
   :fn nil
   :fnargs nil]
  (let [seqin       (fs/fullpath in)
        blast-out   (or (and out (fs/fullpath out))
                        (str seqin ".blast"))
        seq-blastdb (or (and blastdb (fs/fullpath blastdb))
                        default-binary-db)

        blast-path  (get-tool-path :ncbi)
        [blastn
         tblastn]   (map #(str blast-path %) '("blastn" "tblastn"))

        io-args     ["-query" seqin "-db" seq-blastdb "-out" blast-out]
        ctrl-args   ["-strand" (name strand)
                     "-word_size" (str word-size)
                     "-evalue" (str evalue)]
        fmt-args    ["-outfmt" fmt]]

    (assert-tools-exist [blastn tblastn])

    (case pgm
     :blastn
     (runx blastn (concat io-args ctrl-args fmt-args))

     :tblastn
     (let [ctrl-args (vec (drop 2 ctrl-args))]
       (runx tblastn (concat io-args ctrl-args fmt-args))))

    (if-not fn
      blast-out
      (do
        (assert (fn? fn))
        (let [fnargs (if (sequential? fnargs) fnargs (list fnargs))]
          (apply fn blast-out (seq fnargs)))))))


(defn blastn
  "NCBI blastn (from the blast+ toolset).

   IN is the input filespec in fasta format.

   KEYWORD ARGS:
   :out is the output filespec for direct bast output, default `in`.blast
   :blastdb the filespec for the preformatted blast database to use
   :strand (string) plus or minus, defaults to plus
   :word-size blast bp slice size (chunk size for alignment), default 8
   :evalue (real) expectation (E) value threshold for saving hits, default 10
   :fmt NCBI blastn out fmt string, defaults to:
        10 qseqid qstart qend evalue sseqid sstart send
   :fn post processing fn of at least 1 arg for blast output, eg parse-blast
   :fnargs any other args for fn, e.g., output filespec"
  [in & kw-args]
  (apply blast :blastn in kw-args))


(defn tblastn
  "NCBI tblastn (from the blast+ toolset).

   IN is the input filespec in fasta format.

   KEYWORD ARGS:
   :out is the output filespec for direct blast output, default `in`.blast
   :blastdb the filespec for the preformatted blast database to use
   :word-size blast bp slice size (chunk size for alignment), default 8
   :evalue (real) expectation (E) value threshold for saving hits, default 10
   :fmt NCBI blastn out fmt string, defaults to:
        10 qseqid qstart qend evalue sseqid sstart send
   :fn post processing fn of at least 1 arg for blast output, eg parse-blast
   :fnargs any other args for fn, e.g., output filespec"
  [in & kw-args]
  (apply blast :tblastn in kw-args))


(defn blastpgm
  "Determine blast program to use based on sequence alphabet of FASTA-FILE.
   If amino acid alphabet, use tblastn, if nucleotide use blastn.  Returns
   the function for each of these."
  [fasta-file]
  (with-open [r (io/reader fasta-file)]
    (let [l (doall (second (line-seq r)))]
      (if (some #(not (in (first (clojure.string/lower-case %)) "atgcu")) l)
        tblastn
        blastn))))




;;; CD Hit Redundancy removal.
;;;
(defn get-cd-hit-infile [x]
  (if (coll? x)
    (nms-sqs->fasta-file
     (map (fn[[nm rng sq]] [(str nm ":" rng) sq]) x)
     (fs/tempfile "cdhit-" ".fna"))
    (fs/fullpath x)))

(defn cd-hit-est [input outfile & {c :c n :n :or {c 0.90 n 8}}]
  (let [infile (get-cd-hit-infile input)
        outfile (fs/fullpath outfile)
        cdhit-path (get-tool-path :cdhit)
        cdhitcmd (str cdhit-path "cd-hit-est")
        cmdargs ["-i" infile "-o" outfile
                 "-c" (str c) "-n" (str n)]]
    (assert-tools-exist [cdhitcmd])
    (catch-all (runx cdhitcmd cmdargs))
    (when (fs/exists? (str outfile ".clstr"))
      (fs/rm (str outfile ".clstr")))
    (when (fs/exists? (str outfile ".bak.clstr"))
      (fs/rm (str outfile ".bak.clstr")))
    (when (coll? input) (fs/rm infile))
    outfile))




;;; ----------------------------------------------------------------------
;;; CMFinder tools and operations

(defn process-blast-for-cmfinder [acc line]
  "A simple reducer for blastn output where the blast run specified output
   format with the options:
     -outfmt '10 qseqid qstart qend evalue sseqid start send'"
  (let [bs (str/split #"," line)
        qs (take 4 bs)
        ss (drop 4 bs)]
    (conj
     acc
     {:qs (conj (rest qs) (str/replace-re #"^gi\|\d+\|" "" (first qs)))
      :ss (conj (rest ss) (str/replace-re #"^gi\|\d+\|" "" (first ss)))})))




(defn parse-blast [in out]
  "Parse a (NCBI) blastn output for use by cmfinder tools.  The blast run
   specified output format with the options:
     -outfmt '10 qseqid qstart qend evalue sseqid start send <& others>'
   IN is the blast output filespec (path or file obj) and OUT is a filespec
   (path or file obj) of where to place the parsed output.  This file is then
   intended to be use as input to the CANDS program of cmfinder that generates
   per motif-loop data files for input to cmfinder."
  (do-text-to-text
   [in out]
   (let [bs (str/split #"," $line)
         qseq (take 4 bs)
         sseq (take 3 (drop 4 bs))
         [qid qs qe e] qseq
         [sid ss se] sseq]
     (when (and (not (= qid sid)) (< (Float. e) 0.5))
       (printf "%S :\t %S - %S\t %S\n"
               (str/replace-re #"^gi\|\d+\|" "" qid) qs qe e)
       (printf "%S :\t %S - %S\n\n"
               (str/replace-re #"^gi\|\d+\|" "" sid) ss se))))
  out)




(defnk candf [in out
             :max-cand-seq 40
             :min-len-cand 30
             :max-len-cand 100
             :min-stem-loops 1
             :max-stem-loops 1]

  {:pre [(<= min-stem-loops max-stem-loops)
         (<= max-len-cand 100)
         (<= 1 min-len-cand)
         (<= min-len-cand max-len-cand)]}

  (let [seqin        (fs/fullpath in)
        stem-loops   (if (= max-stem-loops min-stem-loops)
                       (str min-stem-loops)
                       (str min-stem-loops "." max-stem-loops))
        seqout       (or out (str seqin ".h" stem-loops ".cand"))

        cmf-path (get-tool-path :cmfinder)
        candf (str cmf-path "candf")

        args (map str ["-c" max-cand-seq "-o" seqout
                       "-M" max-len-cand "-m" min-len-cand
                       "-s" min-stem-loops "-S" max-stem-loops
                       seqin])]
    (assert-tools-exist [candf])
    (runx candf args)))




(defnk cands [in candf-output
             :max-out-motifs 3
             :expected-motif-freq 0.80
             :blast-seq-match nil]
  (let [seqin (fs/fullpath in)
        cmf-path (get-tool-path :cmfinder)
        cands (str cmf-path "cands")
        args (keep #(when % (str %))
                   (flatten ["-n" max-out-motifs "-f" expected-motif-freq
                             (when blast-seq-match ["-m" blast-seq-match])
                             seqin candf-output]))]
    (assert-tools-exist [cands])
    (runx cands args)))




(defn canda [seqin cands-file canda-ofile]
  (let [cmf-path (get-tool-path :cmfinder)
        canda (str cmf-path "canda")]
    (assert-tools-exist [canda])
    (runx canda
          (fs/fullpath cands-file)
          (fs/fullpath seqin)
          (fs/fullpath canda-ofile))))




(defnk cmfinder [seqin canda-file motif-out cm-out
                 :initial-cm nil
                 :candf-output nil
                 :stdout nil]
  (let [cmf-stdout   (and stdout (fs/fullpath stdout))
        initial-cm   (and initial-cm (fs/fullpath initial-cm))
        candf-output (and candf-output (fs/fullpath candf-output))

        cmf-path (get-tool-path :cmfinder)
        cmfinder (str cmf-path "cmfinder")

        args (keep #(when % %)
                   (flatten ["-o" (fs/fullpath motif-out)
                             "-a" (fs/fullpath canda-file)
                             (when initial-cm ["-i" initial-cm])
                             (when candf-output ["-c" candf-output])
                             (fs/fullpath seqin)
                             (fs/fullpath cm-out)
                             (when cmf-stdout [:> cmf-stdout])]))]
    (assert-tools-exist [cmfinder])
    (runx cmfinder args)))




(defnk cmfinder* [in
                  :blastpgm nil
                  :blastdb nil
                  :initial-cm nil
                  :strand :plus
                  :word-size 8
                  :max-cand-seq 40
                  :min-len-cand 30
                  :max-len-cand 100
                  :max-out-motifs 3
                  :expected-motif-freq 0.80
                  :min-stem-loops 1
                  :max-stem-loops 1
                  :del-files true]

  (let [seqin        (fs/fullpath in)

        stem-loops   (if (= max-stem-loops min-stem-loops)
                        (str min-stem-loops)
                        (str min-stem-loops "." max-stem-loops))
        seqout-candf (str seqin ".h" stem-loops ".cand")
        seq-match    (str seqin ".match")]

    ;; CANDF --> produce stem loop candidates (all in seqout-candf)
    (candf seqin seqout-candf
           :max-cand-seq max-cand-seq
           :min-len-cand min-len-cand :max-len-cand max-len-cand
           :min-stem-loops min-stem-loops :max-stem-loops max-stem-loops)

    ;; When blasting, produce blast match file SEQ-MATCH, for input to cands
    (when blastpgm
      (blast  blastpgm seqin (str seqin ".blast")
              :blastdb blastdb
              :strand strand :word-size word-size
              :fn parse-blast :fnargs seq-match))

    ;; CANDS --> produce MAX-OUT-MOTIF candidate files for CANDA/CMFINDER
    (cands seqin seqout-candf
           :max-out-motifs max-out-motifs
           :expected-motif-freq expected-motif-freq
           :blast-seq-match (when blastpgm seq-match))

    ;; Now, for each motif candidate MC, canda MC ... | cmfinder ...
    ;; which yields an alignment and motif (hits) set of output.
    (loop [i (int 0)
           motif-stos []]
      (if (= i max-out-motifs)
        motif-stos
        (let [dot-i (str "." (+ 1 i))
              suffix (str stem-loops dot-i)
              cands-file (str seqout-candf dot-i)]
          (when (not (fs/empty? cands-file))
            (let [canda-ofile (str seqin ".align-sto.h" suffix)
                  motif-ofile (str seqin ".motif-sto.h" suffix)
                  cm-file     (str seqin ".cm.h"    suffix)
                  cmf-stdout  (str seqin ".h" stem-loops ".out" dot-i)]
              (canda seqin cands-file canda-ofile)
              (cmfinder seqin canda-ofile
                        motif-ofile cm-file
                        :initial-cm initial-cm
                        :stdout cmf-stdout)
              (when del-files
                (fs/rm cands-file)
                (fs/rm canda-ofile)
                (fs/rm cm-file)
                (fs/rm cmf-stdout))
              (recur (unchecked-inc i)
                     (conj motif-stos motif-ofile)))))))
    ))




(defn cmalign [cm seqfile outfile
               & {opts :opts par :par :or {opts ["-q" "-1"] par 3}}]
  (let [infernal-path (get-tool-path :infernal)
        cmaligncmd (str infernal-path "cmalign")
        cmfile (fs/fullpath cm)
        seqfile (fs/fullpath seqfile)
        outfile (fs/fullpath outfile)
        mpirun "mpirun"
        cmdargs (conj (vec (concat ["-np" (str par)
                                    cmaligncmd "--mpi" "-o" outfile]
                                   opts))
                      cmfile seqfile)]
    (if (fs/empty? seqfile)
      (raise :type :empty-seq-file :file seqfile)
      (do (assert-tools-exist [cmaligncmd])
          (runx mpirun cmdargs)))
    outfile))

