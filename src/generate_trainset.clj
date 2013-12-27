(ns generate-trainset
  (:require [clojure.string :as str]
            [clojure.contrib.io :as io]
            [edu.bc.fs :as fs])
  (:use [clojure.set :only (map-invert)]
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.tools
        edu.bc.bio.sequtils.snippets-files
        [edu.bc.utils :only (nCk transpose)]
        refold
        ))


(defn split-str-at [s re]
  (map str (re-seq re s) (rest (str/split s re))))

(defn degap-col
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
          sqs (degap-col sqs)))
  ([sqs st]
     (let [dgap (degap-col sqs st)
           st (last dgap)]
       (conj (map (fn [[nm sq] ungap-sq]
                     ;;recombine degapped seqs with the proper names
                     [nm ungap-sq]) 
                   sqs dgap)
             st))))

(defn randnseqs
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

(defn sto->nsubset-stos
  "takes a sto input file and generates n subset sto that contains a
  subset (3-6) of the sequences from the original sto. An aln file is
  first made, then the aln is converted to sto format. Returns a list
  of the file names of the temp stos"

  [insto n]
  (let [[_ seq-lines cons-lines] (join-sto-fasta-lines insto "")
        subset (->> (repeatedly #(-> (randnseqs seq-lines (rand-nth (range 3 7)))
                                     degap-aln))
                    (filter #(apply distinct? %) ) ;elements of subset distinct
                    distinct) ;subsets are distinct
        ]
    (->> (map (fn [ss]
                (let [tmp1 (fs/tempfile)
                      tmp2 (fs/tempfile)]
                  (io/with-out-writer tmp1 (print-aln ss)) ;subset aln
                  (aln->sto tmp1 tmp2) ;subset sto
                  (fs/rm tmp1)
                  (if (valid-seq-sto? tmp2)
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

  [indir outdir file-type]
  (doseq [f (fs/directory-files indir file-type) ;;use trainset/list2 or trainset/neg/list-sto-only
          ]
    (let [c (atom 0)
          seq-lines (read-seqs f)
          fnm (fs/basename f)
          ;;f "RF01693-seed.sto"
          ]
      (if (>= (nCk (count seq-lines) 3) 10)
        (doseq [subset (sto->nsubset-stos f 10)]
          (prn c f)
          (fs/copy subset (str outdir (subs fnm 0 (- (count fnm) 3)) @c ".sto"))
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
  subset sto will have the same structure as the original but the
  gapped columns will be removed."

  [insto n & {:keys [outfile]}]
  (let [[_ seq-lines cons-lines] (join-sto-fasta-lines insto "")
        st (->> cons-lines
                (filter #(.startsWith (first %) "#=GC SS_cons") )
                (map #(-> % second last)  )
                first
                change-parens)
        [outst & outseqs] (-> (randnseqs seq-lines n)
                              (degap-aln st))]
    (if outfile
      (io/with-out-writer outfile (print-sto outseqs outst))
      (print-sto outseqs outst))
    ))

(defn entry-parts->entry
  "changes a entry vector back into the string"

  [entry-vector]
  (let [[nm [s e] strand] entry-vector]
    (str/join "/" [nm (str s "-" e) strand])))

(defn get-flank
  "Given an entry it will find the flanking 5' and 3' region of the
  same size as the entry. Flanking regions are returned as a map with
  the keys :5prime and :3prime. "

  [entry]
  (let [[_ [s e] strand] (entry-parts entry)
        delta (- e s)
        lentry (gen-name-seq entry :rdelta (- delta) :ldelta delta)
        rentry (gen-name-seq entry :rdelta delta :ldelta (- delta))
        flank (if (pos? (Integer/parseInt strand)) ;make strand an int
                {:5prime lentry :3prime rentry}
                {:5prime rentry :3prime lentry})]
    flank))

(defn trainset3-pos
  "Takes the Rfam stos and produces stos with NC assession
  numbers. Any seq that doesn't have a NC equivalent is left out of
  the new sto" []
  (map embl-to-nc (fs/directory-files "/home/peis/bin/gaisr/trainset3/pos" ".sto")))

(def cur-struct-negs
  '("RF00514-seed-NC-3prime.sto" "RF01055-seed-NC-3prime.sto" "RF01826-seed-NC-3prime.sto" "RF00558-seed-NC-3prime.sto" "RF01070-seed-NC-3prime.sto" "RF01510-seed-NC-3prime.sto" "RF01693-seed-NC-3prime.sto" "RF00559-seed-NC-3prime.sto" "RF01769-seed-NC-3prime.sto" "RF01727-seed-NC-3prime.sto" "RF01385-seed-NC-3prime.sto" "RF01692-seed-NC-3prime.sto" "RF01402-seed-NC-3prime.sto" "RF01482-seed-NC-3prime.sto" "RF01767-seed-NC-3prime.sto" "RF01694-seed-NC-3prime.sto" "RF00506-seed-NC-3prime.sto" "RF00114-seed-NC-3prime.sto" "RF01402-seed-NC-5prime.sto" "RF01482-seed-NC-5prime.sto" "RF01769-seed-NC-5prime.sto" "RF01826-seed-NC-5prime.sto" "RF01727-seed-NC-5prime.sto" "RF01694-seed-NC-5prime.sto" "RF01693-seed-NC-5prime.sto" "RF01510-seed-NC-5prime.sto" "RF01385-seed-NC-5prime.sto"))

(defn trainset3-negs
  "generate trainset negative sto's. Produces the stos in current
  folder. Usually just copy to a new folder for negatives"

  []
  (let [fdir "/home/peis/bin/gaisr/trainset3/pos"
        files (fs/re-directory-files fdir "-NC.sto")
        outfn (fn [k f outfile-ext] ;outfile-ext is really the outfile
                                   ;name part to appear
                (let [outfile (fs/replace-type f outfile-ext)
                      flank-entries (map get-flank (read-seqs f :info :name))]
                  (-> (nms-sqs->fasta-file (-> (map k flank-entries);k=:5prime or :3prime
                                               map-invert 
                                               map-invert ;removes duplicates
                                               vec)
                                           (fs/replace-type outfile ".fna")) ;fasta out
                      fasta->aln
                      (aln->sto outfile))))] ;out sto file
    (doseq [f files]
      (outfn :5prime f "-5prime.sto")
      (outfn :3prime f "-3prime.sto"))))

(defn trainset3-shuffled []
  (let [fdir "/home/peis/bin/gaisr/trainset3/"
        files (fs/re-directory-files (str fdir "pos") #"-NC.sto$")]
    (map (fn [f]
           (let [outname (fs/join fdir
                                  "shuffled"
                                  (-> (fs/replace-type f ".shuffle.sto")
                                      fs/basename))]
             (sto->randsto f outname 10)))
         files)))

(defn check-blastout-hits
  "checks blastout hits to identify hits that are in the blastout file
  that appear to be 'valid' hits but are not due to incorrect length. "

  [stoin]
  (let [blastout (fs/replace-type stoin ".fna.blast")
        embl-aln-pairs (->> stoin (#(read-seqs % :info :both)) (into {}))
        ;;compare the ncbi seq and the inseq from alignment and the
        ;;seq lengths to ensure they are the same
        embl-ncbi-seq-len (fn [[embl-name ncbi-hits]] 
                            (map (fn [entry]
                                   (let [[nm [s e] strand] (entry-parts entry)
                                         nc-seq (second (gen-name-seq entry))
                                         rfam-seq (-> (embl-aln-pairs embl-name)
                                                       (str/replace #"\." "" ))]
                                     [nm
                                      (= (- e s -1) (count rfam-seq))
                                      (= nc-seq rfam-seq)]))
                                 ncbi-hits))
        good-mappings (->> (map #(vector (first %) (embl-ncbi-seq-len %))
                                (get-embl-blast-candidates blastout))
                           (map (fn [[embl-name entries]]
                                  (some #(and (-> % second true?)
                                              (->> % last true?))
                                        entries)))
                           (group-by true? ))]
    (/ (count (good-mappings false))
       (+ (count (good-mappings true))
          (count (good-mappings false))))
    ))

(defn check-sto-entry-seqs
  "checks all entries in the sto to make sure the lengths are all correct."
  
  [stoin]
  (every? true? (for [[entry s] (read-seqs stoin :info :both)
                      :let [length (fn [[s e]] (- e s -1))
                            coord (second (entry-parts entry))]]
                  (= (length coord) (count (str/replace s #"\." ""))))))
