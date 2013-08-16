(ns edu.bc.bio.sequtils.snippets-files
  (:require [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.java.shell :as shell]
            [clojure.set :as set]
            [edu.bc.fs :as fs])
  (:use [clojure.contrib.pprint :only (cl-format)]
        edu.bc.utils.fold-ops
        edu.bc.bio.sequtils.files
        [edu.bc.utils.probs-stats :only (probs)]
        [refold :only (remove-gaps)]
        [slingshot.slingshot :only [throw+]]))

(defn change-parens
  "Change the stucture line to use (, ), and  isntead of <,>,-,and :"
  
  [struct]
  (->> struct
       (str/replace-re #"\<" "(") ;open
       (str/replace-re #"\>" ")" ) ;close
       (str/replace-re #"\:|\-" "." ))) ;replace gaps

(defn read-sto
  "read stockholm file creates a map where the key=name
  val=sequence. the keys are :seqs :cons :file :cov"
  
  [f & {:keys [info] 
        :or {info :data}}]
  (let [[gc-lines seq-lines cons-lines] (join-sto-fasta-lines f "")
        cov (first (map #(last (str/split #"\s+" %))
                        (filter #(.startsWith % "#=GC cov_SS_cons") gc-lines)))
        cl (map #(last (second %))
                (filter #(.startsWith (first %) "#=GC SS_cons") cons-lines))
        sl (reduce (fn [v [nm [_ sq]]]
                     (let [sq (str/replace-re #"T" "U" (.toUpperCase sq))]
                       (case info
                         :both (conj v [nm sq])
                         :name (conj v nm)
                         :data (conj v sq))))
                   [] seq-lines)]
    {:seqs sl :cons cl :file f :cov cov}))

(defn read-aln
  "read clustalW file and makes it unblocked"

  [aln & {:keys [info]
          :or {info :data}}]
  (let [seq-lines (second (join-sto-fasta-lines aln ""))
        sl (reduce (fn [v [nm [_ sq]]]
                     (let [sq (str/replace-re #"T" "U" (.toUpperCase sq))]
                       (case info
                        :both (conj v [nm sq])
                        :name (conj v nm)
                        :data (conj v sq))))
                   [] seq-lines)]
    (-> sl rest butlast vec)))

(defn- valid-seq-struct
  "Checks the sqs (list of sequences) to make sure that all sequences
   can form part of the consensus structure. This is useful when
   ensuring that the shuffled stos will form valid
   structures. Occassionally, sequences will not fold leaving an empty
   structure and causing other functions to fail. Returns true if all
   sequences can fold into part of the cons structure."
  
  [sqs cons]
  (let [valid? (fn [st] (pos? (count (struct->matrix st))))]
    (every? true? (map (fn [[_ s]]
                         (let [[_ st] (remove-gaps s cons)]
                           (valid? st)))
                       sqs))))

(defn valid-seq-sto?
  "Checks the sto to make sure that all sequences in the file can form
   part of the consensus structure. This is useful when ensuring that
   the shuffled stos will form valid structures. Occassionally,
   sequences in the sto will not fold into the consensus structure and
   cause other functions to fail. Returns true if all sequences can
   fold into part of the cons structure."
  
  [sto]
  (let [{sqs :seqs cons :cons} (read-sto sto :info :both)
        cons (change-parens (first cons))]
    (valid-seq-struct sqs cons)))

(defn print-aln
  "prints out the seq-lines in clustalW format"
  
  [seq-lines]
  (println "CLUSTAL W (1.83) multiple sequence alignment\n")
  (doseq [[nm sq] seq-lines]
    (cl-format true "~A~40T~A~%" nm sq)))

;;;(ns-unmap 'edu.bc.bio.sequtils.snippets-files 'fasta->aln)
(defmulti

  ^{:doc "align a fasta file using either :clustalw
  or :centroid_align. Produces a file aln-out. If no outfile name is
  specified then the current file name is used with the extension
  changed to '.aln'.

  TODO: should add muscle as an aligner."
    :arglists '([alignprogram fasta-in] [alignprogram fasta-in aln-out])}

  fasta->aln (fn [& args]
               (first args)))

(defmethod fasta->aln :clustalw [_ fasta-in aln-out]
  (let [call (shell/sh "clustalw"
                       (str "-infile=" fasta-in)
                       (str "-outfile=" aln-out)
                       "-quiet")]
    (if-not (empty? (call :err))
      (throw+ {:error (call :err)})
      aln-out)))

(defmethod fasta->aln :centroid_align [_ fasta-in aln-out]
  (let [call (shell/sh "centroid_align"
                       "-f" "clustalw"
                       fasta-in)]
    (if-not (empty? (call :err))
      (throw+ {:error (call :err)})
      (io/with-out-writer aln-out (-> call :out print)))
    aln-out))

(defmethod fasta->aln :default [alignprogram fasta-in]
  (fasta->aln alignprogram fasta-in (fs/replace-type fasta-in ".aln")))

(comment
  (defn fasta->aln
   "aligns sequences in a fasta file using clustalW. Alignment
  automatically printed by clustalW to same folder location containing
  .aln extension.

  TODO: add multimethod for using muscle and centroid_align"

   ([fasta-in]
      (fasta->aln fasta-in (fs/replace-type fasta-in ".aln")))
  
   ([fasta-in aln-out]
      (let [call (shell/sh "clustalw"
                           (str "-infile=" fasta-in)
                           (str "-outfile=" aln-out)
                           "-quiet")]
        (if-not (empty? (call :err))
          (throw+ {:error (call :err)})
          aln-out)))))

(defn- print-sto
  "takes sequence lines and a structure line and writes it into a sto
  format file. the seq-lines needs to be a collection of [name
  sequence] pairs. structure is a string. Simply prints out to the
  repl."

  ([seq-lines structure]
     (println "# STOCKHOLM 1.0\n")
     (doseq [sq seq-lines]
       (let [[nm sq] (if (vector? sq)
                       sq
                       (str/split #"\s+" sq))]
         (cl-format true "~A~40T~A~%" nm (str/replace-re #"\-" "." sq))))
     (cl-format true "~A~40T~A~%" "#=GC SS_cons" structure)
     (println "//")))

;;;(ns-unmap 'edu.bc.bio.sequtils.snippets-files 'aln->sto)
(defmulti aln->sto (fn [in-aln out-sto & args]
                     (let [args (first args)]
                       [(type (args :st)) (args :foldmethod)])))

(defmethod aln->sto [String ::RNAalifold] [in-aln out-sto args]
  (let [sq (read-aln in-aln :info :both)]
    (io/with-out-writer out-sto
      (print-sto sq (args :st)))
    out-sto))

(defmethod aln->sto [nil ::RNAalifold] [in-aln out-sto args]
  (let [st (fold-aln in-aln)
        sq (read-aln in-aln :info :both)]
    (io/with-out-writer out-sto
      (print-sto sq st))
    out-sto))

(defmethod aln->sto [nil ::centroid_alifold] [in-aln out-sto args]
  (let [st (fold-aln in-aln {:foldmethod ::centroid_alifold})
        sq (read-aln in-aln :info :both)]
    (io/with-out-writer out-sto
      (print-sto sq st))
    out-sto))

(defmethod aln->sto :default [in-aln out-sto]
  (aln->sto in-aln out-sto {:foldmethod ::RNAalifold}))

(comment 
 (defn aln->sto
   "takes an alignment in Clustal W format and produces a sto file by
   using RNAalifold to determine the structure and then making it into
   a sto file adding header and a consensus line"

   [in-aln out-sto & {:keys [fold-alg st]
                      :or {fold-alg "RNAalifold"}}]
   (cond
    (identity st)                       ;structure provided
    (let [sq (read-aln in-aln :info :both)]
      (io/with-out-writer out-sto
        (print-sto sq st))
      out-sto)

    (= fold-alg "RNAalifold")           ;structure from RNAalifold
    (let [st (fold-aln in-aln)
          sq (read-aln in-aln :info :both)]
      (io/with-out-writer out-sto
        (print-sto sq st))
      out-sto)                          ;return out sto filename

    ;;else use cmfinder
    #_(shell/sh "perl" "/home/kitia/bin/gaisr/src/mod_cmfinder.pl" in-aln out-sto))))

(defn sto->randsto
  "takes a sto input file and generates a random sto to specified
  output. Program uses SISSIz to make the random aln. Then the aln is
  converted to sto format by adding a structure line. Structure is
  found by using RNAalifold. Returns the file name of the new sto"

  [insto outsto]
  (let [fname (str/take (- (count insto) 3) insto)
        inaln (sto->aln insto (str fname "aln"));convert in-sto to aln
        outaln (str fname "out")]
    (io/with-out-writer outaln
      (print ((shell/sh "SISSIz" "-s" inaln) :out))) ;create random aln
    (aln->sto outaln outsto)
    ;;remove temp files
    (fs/rm inaln)
    (fs/rm outaln)
    outsto))

(defn sto->fasta
  "converts a sto file to a fasta file. removes the gaps from the
   sequences. This should be used in order to make a fasta file for
   input into CMfinder"
  
  ([sto]
     (let [lines (second (join-sto-fasta-lines sto ""))]
       (doseq [[nm [_ sq]] lines]
         (println (str ">" nm))
         (println (str/replace-re #"\." "" sq)))))
  
  ([sto outfasta]
     (io/with-out-writer outfasta
      (sto->fasta sto))))



(defn profile
  "Creates a profile object of the sto which is read in. The keys
   include sequences, structure, base-frequency of columns, background
   base frequency, which columns are paired, covariation, filename and
   length."
  
  [sto]
  (let [p (read-sto sto)
        struct (map #(change-parens %) (p :cons))
        s (p :seqs)
        freqs (partition 2 (interleave (range (count (first s)))
                                       (map frequencies
                                            (apply map vector (map #(rest (str/split #"" %)) s)))))
        fract-freqs (sort-by key (reduce (fn [l [n freq-map]]
                                           (assoc l n (probs freq-map)))
                                         {} freqs))
        q (probs 1 (rest (str/split #"" (apply str (p :seqs)))))
        pairs (map #(refold/make_pair_table %) struct)]
  {:seqs s
   :structure struct
   :fract fract-freqs
   :background q
   :pairs pairs
   :cov (p :cov)
   :filename (p :file)
   :length  (count (first s))}))

(defn read-clj
  "Reads a clj data structure"
  
  [f]
  (->> (slurp f) read-string ))
