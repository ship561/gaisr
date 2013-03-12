(ns generate-trainset
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])
  (:use edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.snippets-files
        [edu.bc.utils :only (nCk transpose)]
        edu.bc.bio.sequtils.snippets-files
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))

(defn randseqs
  "takes seq-lines and the number of random sequences to draw from it"

  [seq-lines n]
  (let [rand-sqs (->> (repeatedly 1000 #(shuffle seq-lines))
                      last
                      (take n))
        remove-gap-col (fn [cols]
                         (->> (map #(-> % second second) cols) ;gets
                                        ;seqs
                              transpose
                              (remove #(empty? (str/replace-re #"\.|\-" "" %)))
                              transpose))]
    (->> (map (fn [[nm [uid sq]] ungap-sq]
                [nm [uid ungap-sq]]) ;recombine degapped seqs
                                        ;with the proper names
              rand-sqs
              (remove-gap-col rand-sqs))
         (sort-by #(-> % second first) ) ;preserve order of seqs from original
         (map (fn [[nm [_ sq]]] [nm sq]) );only return name and seq
         )))

(defn sto->subset-sto
  "takes a sto input file and generates a random sto that contains a
  subset (n) of the sequences in the original sto. An aln file is
  first made, then the aln is converted to sto format. Returns a list
  of the file names of the temp stos"

  [insto n]
  (let [[_ seq-lines _] (join-sto-fasta-lines insto "")
        subset (->> (repeatedly #(randseqs seq-lines (+ (rand-int 4) 3)))
                    (filter #(apply distinct? %) ) ;elements of subset distinct
                    distinct) ;subsets are distinct
        ]
    (->> (map (fn [ss]
                (let [tmp1 (fs/tempfile)
                      tmp2 (fs/tempfile)]
                  (io/with-out-writer tmp1 (toaln ss)) ;subset aln
                  (aln->sto tmp1 tmp2) ;subset sto
                  (fs/rm tmp1)
                  (if (valid-seq-sto tmp2)
                    tmp2               ;valid sto
                    (do (fs/rm tmp2)
                        :remove))))
              subset)
         (remove #(= :remove %) ) 
         (take n)) ;take first n valid stos
    ))

(defn make-training-set
  "make the training set from existing stos. The training set is
  composed of stos which are a subset of the positive or negative
  stos. The positive stos are from Rfam and the negative stos are
  generated using SISSIz to preserve mono/dinucleotide composition as
  well as column base conservation and gap conservation. The indir is
  the directory of stos and the outdir is where the training stos are
  put.

  Usually use trainset/neg or trainset/ (for pos)."

  [indir outdir]
  (doseq [f (->> indir
                 (fs/listdir )
                 (filter #(re-find #".sto$" %) )) ;;use trainset/list2 or trainset/neg/list-sto-only
          ]
    (let [c (atom 0)
          seq-lines (read-seqs (str indir f))
          ;;f "RF01693-seed.sto"
          ]
      (if (>= (nCk (count seq-lines) 3) 10)
        (doseq [subset (sto->subset-sto (str indir f) 10)]
          (prn c f)
          (fs/copy subset (str outdir (subs f 0 (- (count f) 3)) @c ".sto"))
          (swap! c inc)
          (fs/rm subset))
        (println f "did not work. too few sequences")
        ))))

(defn foo
  "don't remember what this does yet."
  
  []
  (let [d "/home/kitia/bin/gaisr/trainset2/"
        ]
    (doseq [stos (fs/listdir d)
            s stos]
      (let [fasta (str (subs s 0 (- (count s) 3)) "fasta")]
        (sto->fasta (str d s) (str d fasta))))))




