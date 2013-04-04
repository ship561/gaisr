(ns generate-trainset
  (:require [clojure.string :as str]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])
  (:use edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.snippets-files
        [edu.bc.utils :only (nCk transpose)]
        edu.bc.bio.sequtils.snippets-files
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        refold
        ))

(defn degap
  "Takes sqs which is a coll of seqs [name sequence] from an alignment
  with/without a structure st. Removes all gapped colls and if a
  structure is present, adjusts the structure to reflect the removed
  gaps. Returns alignment as a vector of strings. Last element is the
  structure."

  ([sqs]
     (->> (map second sqs) ;gets seqs
          transpose
          (remove #(empty? (str/replace % #"\.|\-" "")))
          transpose))

  ([sqs st]
     (let [;sqs ["A.AA...UUU" "AAA....UUU"]
           sqs (mapv second sqs)
           ;st (first st)
           ;;st "((......))"
           table (when st (transient (make_pair_table st)))
           pairs #{"AU" "UA" "CG" "GC" "GU" "UG"}
           len (count st)]
       (->> (reduce (fn [x alncoll]
                      (let [l (count alncoll)
                            seqs (butlast alncoll)
                            i (last alncoll)]
                        (conj x
                              (cond
                               (every? #(= \. %) seqs) ;all gaps
                               (repeat (count seqs) \x)
                               
                               ;;seqs 
                               (and (table i)
                                    (->> (map (fn [c inseq]
                                                (->> (str c (.charAt inseq (table i)))
                                                     (contains? pairs )))
                                              (butlast seqs) sqs)
                                         (not-any? true? )))
                               (do (dissoc! table (table i) i)
                                   (repeat (count seqs) \x))
                               
                               (not (table i)) ;not bp location turns into gap
                               (conj (vec (butlast seqs)) \.)
                               
                               :else seqs))))
                    [] (apply map vector (conj sqs st (range len))))
            transpose
            (map (fn [invec] (remove #(= \x %) invec))) ;removes x'd gaps
            (map str/join )))))

(defn degap-aln
  
  ([sqs]
     (map (fn [[nm sq] ungap-sq]
            ;;recombine degapped seqs with the proper names
            [nm ungap-sq]) 
          sqs (degap sqs)))
  ([sqs st]
     (let [dgap (degap sqs st)
           st (last dgap)]
       (conj (map (fn [[nm sq] ungap-sq]
                     ;;recombine degapped seqs with the proper names
                     [nm ungap-sq]) 
                   sqs dgap)
             st))))

(defn randseqs
  "takes n seq-lines and the number of random sequences to draw from it"

  [seq-lines n]
  (let [rand-sqs (->> (repeatedly 1000 #(shuffle seq-lines))
                      last
                      (take n))
        ]
    (->> rand-sqs
         (sort-by #(-> % second first) ) ;preserve order of seqs from original
         (map (fn [[nm [_ sq]]] [nm sq]) );only return name and seq
         )))

(defn sto->subset-sto
  "takes a sto input file and generates n subset sto that contains a
  subset (3-6) of the sequences from the original sto. An aln file is
  first made, then the aln is converted to sto format. Returns a list
  of the file names of the temp stos"

  [insto n]
  (let [[_ seq-lines cons-lines] (join-sto-fasta-lines insto "")
        subset (->> (repeatedly #(-> (randseqs seq-lines (+ (rand-int 4) 3))
                                     degap-aln))
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





(defn make-subset-sto
  "Takes a sto and randomly chooses n sequences to make a new sto. The
  subset sto will have the same structure as the original minus the
  gapped columns will be removed."

  [insto n & {:keys [outfile]}]
  (let [[_ seq-lines cons-lines] (join-sto-fasta-lines insto "")
        st (->> cons-lines
                (filter #(.startsWith (first %) "#=GC SS_cons") )
                (map #(-> % second last)  )
                first
                change-parens)
        [outst & outseqs] (-> (randseqs seq-lines n)
                              (degap-aln st))]
    (if outfile
      (io/with-out-writer outfile (print-sto outseqs outst))
      (print-sto outseqs outst))
    ))
