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
                         [nm [uid ungap-sq]])
                       rand-sqs
                       (remove-gap-col rand-sqs))
         (sort-by #(-> % second first) )
         (map (fn [[nm [_ sq]]] [nm sq]) )
         )))

(defn make-training-set []
 ;;code looks through list of stos and picks random sequences out of
 ;;the sto
  (doseq [f (->> (str (fs/homedir) "/bin/gaisr/trainset/neg")
                 (fs/listdir )
                 (filter #(re-find #".sto$" %) )) ;;use trainset/list2 or trainset/neg/list-sto-only
          ]
    (let [c (atom 0)
          fdir "/home/kitia/bin/gaisr/trainset/neg/" ;;use either trainset/ or trainset/neg
          odir "/home/kitia/bin/gaisr/trainset2/neg/" ;;use either trainset2/pos or trainset2/neg
          seq-lines (read-seqs (str fdir f))
          ;;f "RF01693-seed.sto"
          ]
      (if (>= (nCk (count seq-lines) 3) 10)
        (doseq [subset (sto->subset-sto (str fdir f) 10)]
          (prn c f)
          (fs/copy subset (str odir (subs f 0 (- (count f) 3)) @c ".sto"))
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
