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

(defn- remaining-files [pred]
  (filter pred (fs/listdir (str homedir "/bin/gaisr/trainset2/pos/"))))

(defn create-inv-seqs
  "Generates inverse folded seqs using inverse-fold. If an outfile
   exists, then it will read it in and then add to the existing list
   of seqs. Takes a sequence name(nm), structure (st), n inverse seqs
   to make, and outfile. Returns the list of sequences."

  [nm st n outfile]
  
  (let [cur-seqs (if (fs/exists? outfile)
                   (-> outfile slurp read-string)
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
  
  ([nseqs timeout & {:keys [units ncore]
                     :or {units :s ncore 1}}]
     (let [fdir (str homedir "/bin/gaisr/trainset2/pos/")
           ofile (str homedir "/bin/gaisr/robustness/subopt-robustness-test2.clj")
           done-files (->> (slurp "/home/peis/bin/gaisr/robustness/subopt-robustness-test2.clj") 
                         read-string 
                         (into {}))
           pred (fn [x] (let [outfile (str fdir (str/butlast 3 x) "inv.clj")]
                         (and (re-find #"\.7\.sto" x) ;subset of data
                              (not (contains? done-files (keyword x)))
                              (let [foo (->> (slurp outfile) read-string (into {}))
                                    bar (map count (vals foo))] 
                                (or (< (count (keys foo)) 3) ;correct #seqs
                                    (some #(< % 100) bar)))))) ;correct #invfolds
           diff (take 100 (remaining-files pred))
           timeout-ms (case units
                            :ms timeout
                            :s (* timeout 1000)
                            :min (* timeout 1000 60)
                            :hr (* timeout 1000 60 60))]
       (pxmap (fn [insto]
                (prn "working on file" insto)
                (create-inv-sto (str fdir insto) nseqs timeout-ms))
              ncore
              diff)
       #_(pxmap (fn [instos]
                (doall
                 (for [insto instos
                       :let [outfile (str fdir (str/butlast 3 insto) "inv.clj")]
                       :when (or (not (fs/exists? outfile))
                                 (< (fs/size outfile) 8000))]
                   (do (prn "working on file" insto)
                       (create-inv-sto (str fdir insto) nseqs timeout-ms)))))
              ncore
              diff)
       #_(doall
        (for [instos diff
              insto instos
              :let [outfile (str fdir (str/butlast 3 insto) "inv.clj")]
              :when (or (not (fs/exists? outfile))
                        (< (fs/size outfile) 8000))]
          insto)))))

(defn -main [& args]
  (doall (driver-create-inv 100 10 :units :hr :ncore 5)))

