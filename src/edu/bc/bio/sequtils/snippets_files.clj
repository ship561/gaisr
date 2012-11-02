(ns edu.bc.bio.sequtils.snippets-files
  (:require [clojure.contrib.io :as io]
            [clojure.contrib.string :as str]
            [clojure.java.shell :as shell]
            [edu.bc.fs :as fs])
  (:use edu.bc.bio.sequtils.files
        [clojure.contrib.pprint
         :only (cl-format compile-format)]))

(defn aln->sto
  "takes an alignment in Clustal W format and produces a sto file by using RNAalifold to determine
   the structure and then making it into a sto file adding header and a consensus line"

  [in-aln out-sto & {fold_alg :fold_alg :or {fold_alg "RNAalifold" }}]
  (if (= fold_alg "RNAalifold")
    (let [st (->> ((shell/sh "RNAalifold"
                             "-P" "/home/kitia/Desktop/ViennaRNA-2.0.0/rna_andronescu2007.par"
                             "-r" "--noPS" in-aln) :out)
                  (str/split-lines)
                  second
                  (str/split #" ")
                  first)
          sq (rest (second (join-sto-fasta-lines in-aln "")))]
      (io/with-out-writer out-sto
        (println "# STOCKHOLM 1.0\n")
        (doseq [[n [_ s]] sq]
          (cl-format true "~A~40T~A~%" n (str/replace-re #"\-" "." s)))
        (cl-format true "~A~40T~A~%" "#=GC SS_cons" st)
        (println "//"))
      out-sto) ;return out sto filename
    ;;else use cmfinder
    #_(shell/sh "perl" "/home/kitia/bin/gaisr/src/mod_cmfinder.pl" in-aln out-sto)))



(defn sto->randsto
  "takes a sto input file and generates a random sto to specified
  output. Program uses SISSIz to make the random aln. Then the aln is
  converted to sto format. Returns the file name of the new sto"

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

(defn read-sto
  "read stockholm file creates a map where the key=name val=sequence"
  
  [f & {with-names :with-names
                     :or {with-names false}}]
  (let [[gc-lines seq-lines cons-lines] (join-sto-fasta-lines f "")
        cov (first (map #(last (str/split #"\s+" %))
                (filter #(.startsWith % "#=GC cov_SS_cons") gc-lines)))
        cl (map #(last (second %))
                (filter #(.startsWith (first %) "#=GC SS_cons") cons-lines))
        sl (reduce (fn [v [nm [_ sq]]]
                     (let [sq (str/replace-re #"T" "U" (.toUpperCase sq))]
                       (if-not with-names
                         (conj v sq)
                         (conj v [nm sq]))))
                [] seq-lines)]
    (assoc {} :seqs sl :cons cl :file f :cov cov)))

(defn change-parens
  "Change the stucture line to use (, ), and  isntead of <,>,-,and :"
  
  [struct]
  (->> struct
       (str/replace-re #"\<" "(") ;open
       (str/replace-re #"\>" ")" ) ;close
       (str/replace-re #"\:|\-" "." ))) ;replace gaps
