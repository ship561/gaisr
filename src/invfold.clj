(ns invfold
  (:gen-class)
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str])
  (:use [edu.bc.utils :only (pxmap)]
        edu.bc.utils.fold-ops
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens)]
        refold))

(def ^{:private true} homedir (fs/homedir))

(def ^{:private true} fdir (str homedir "/bin/gaisr/trainset2/pos/"))

(def ^{:private true} todo-files
  (filter #(re-find #"\.7\.sto" %) (fs/listdir fdir))) ;subset of data

(def ^{:private true} done-files
  (let [ofile (str homedir "/bin/gaisr/robustness/subopt-robustness-test2.clj")]
    (->> (read-clj ofile) (into {})  keys)))

#_(defn- remaining-files [outfile]
  (let [ofile outfile ;storage location
        fdir (str homedir "/bin/gaisr/trainset2/pos/")
        done-files (when (fs/exists? ofile) (->> (read-string (slurp ofile)) ;read existing data
                                                 (into {})))]
    (->> (filter #(and (re-find #"\.7\.sto" %) ;subset of data
                       (not (contains? done-files (keyword %)))) ;remove done files
                 (fs/listdir fdir))
         (partition-all 2 ) ;group into manageable
                                        ;chuncks
         )))

(defn- read-clj
  "Reads a clj data structure"
  
  [f]
  (->> (slurp f) read-string ))

(defn- remaining-files
  "Return a list of files to loop over. files todo and files-ignore
   are lists of file names. They need to intersect tof files to be
   removed. Also takes a pred which are usually more specific features
   of the current run."
  
  [pred files-todo files-ignore]
  (filter #(and pred
                (not (contains? (set files-ignore) %)))
          files-todo))

(defn create-inv-seqs
  "Generates inverse folded seqs using inverse-fold. If an outfile
   exists, then it will read it in and then add to the existing list
   of seqs. Takes a sequence name(nm), structure (st), n inverse seqs
   to make, and outfile. Returns the list of sequences."

  [nm st n outfile]
  (let [cur-seqs (if (fs/exists? outfile)
                   (read-clj outfile)
                   {nm []})
        cur-n (count (cur-seqs nm))
        inv-seqs (distinct (lazy-cat (cur-seqs nm) (inverse-fold st n :perfect? false)))]
    (if (>= cur-n n)
      (cur-seqs nm)
      (do (doall inv-seqs)
          (io/with-out-writer outfile
            (prn (assoc cur-seqs nm (vec inv-seqs))))
          (take n inv-seqs)))))

(defn create-inv-sto
  "generates inverse sequences for a sto by calling the
   create-inv-seqs function. Will timeout after timeout-min is
   reached. Input a sto, the number of inverse seqs to generate for
   each seq in the sto and the timeout in minutes. Returns a vector
   [sto-name status] at to indicate success."

  [insto n timeout-ms]
  (let [outfile (str (str/butlast 3 insto) "inv.clj")
        {inseqs :seqs cons :cons} (read-sto insto :with-names true)
        cons (change-parens (first cons))
        f (fn []
            (doall
             (for [[nm s] inseqs]
               (let [[_ st] (remove-gaps s cons)
                     x (create-inv-seqs nm st n outfile)]
                 x)
               )))
        fc-g (future-call f)]
    (if-let [v (deref fc-g timeout-ms false)]
      [insto :done]
      (do (future-cancel fc-g) [insto :cancelled]))))

(defn driver-create-inv
  "drives the create-inv-sto function by feeding it all the stos of
   interest - mainly the *.7.sto. creates nseqs inverse sequences."
  
  ([todo-files done-files nseqs timeout & {:keys [units ncore]
                                           :or {units :s ncore 1}}]
     (let [;;also filter stos that need to be done because they lack
           ;;the correct number of inverse-seqs
           pred (fn [x] (let [outfile (str fdir (str/butlast 3 x) "inv.clj")
                             invseq (->> (read-clj outfile) (into {}))
                             totalinvseq (map count (vals invseq))] 
                         (or (< (count (keys invseq)) 3) ;correct #seqs
                             (some #(< % 100) totalinvseq)))) ;correct #invfolds
           diff (take 100 (remaining-files pred (map keyword todo-files) done-files))
           timeout-ms (case units
                        :ms timeout
                        :s (* timeout 1000)
                        :min (* timeout 1000 60)
                        :hr (* timeout 1000 60 60))]
       (pxmap (fn [insto]
                (prn "working on file" insto)
                (create-inv-sto (str fdir insto) nseqs timeout-ms))
              ncore
              diff))))


(defn -main [& args]
  (let [[todo done] args
        todo (or todo todo-files)
        done (or done done-files)]
    (doall (driver-create-inv todo done 100 10 :units :hr :ncore 5))))

