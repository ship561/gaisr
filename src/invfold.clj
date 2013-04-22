(ns invfold
  (:gen-class)
  (:require [edu.bc.fs :as fs]
            [clojure.contrib.io :as io]
            [clojure.contrib.string :as str])
  (:use [edu.bc.utils :only (pxmap)]
        edu.bc.utils.fold-ops
        [edu.bc.bio.sequtils.snippets-files
         :only (read-sto change-parens read-clj)]
        refold
        [clojure.tools.cli :only [cli]]))

(def ^{:private true} homedir (fs/homedir))

(def ^{:private true} fdir (str homedir "/bin/gaisr/trainset2/pos/"))

(def ^{:private true} todo-files
  (filter #(re-find #"\.7\.sto" %) (fs/listdir fdir))) ;subset of data

(def ^{:private true} done-files
  (let [ofile (str homedir "/bin/gaisr/robustness/subopt-robustness-test2.clj")]
    (if (fs/exists? ofile)
      (->> (read-clj ofile) (into {})  keys (map str/as-str))
      [])))

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

(defn remaining-files
  "Return a list of files to loop over. files todo and files-ignore
   are lists of file names. They need to intersect tof files to be
   removed. Also takes a pred which are usually more specific features
   of the current run."
  
  [pred files-todo files-ignore]
  (remove #(or (pred %)
               (contains? (set files-ignore) %))
          files-todo))

(defn enough-inv-seq?
  [f n]
  (let [invfile (str (str/butlast 3 f) "inv.clj")]
    (if (fs/exists? invfile)
      (let [invseqs (->> (read-clj invfile) (into {}))
            totalinvseqs (map count (vals invseqs))
            cnt-seqs (-> (read-sto f) :seqs count)] 
        (and (>= (count (keys invseqs)) cnt-seqs)      ;correct #seqs
             (every? #(>= % n) totalinvseqs))) ;correct #invfolds
      false)))

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
        inv-seqs (distinct (lazy-cat (cur-seqs nm)
                                     (inverse-fold st n :perfect? false)))]
    (when (< cur-n n)
      (do (doall inv-seqs)
          (io/with-out-writer outfile
            (prn (assoc cur-seqs nm (vec inv-seqs))))))
    (take n inv-seqs)))

(defn create-inv-sto
  "generates inverse sequences for a sto by calling the
   create-inv-seqs function. Will timeout after timeout-ms is
   reached. Input a sto, the number of inverse seqs to generate for
   each seq in the sto and the timeout in milliseconds. Returns a vector
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
           diff (remaining-files #(enough-inv-seq? (str fdir %) nseqs) todo-files done-files)
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


(defn -main
  "args are given as 1 string on the command line with respect to the
   flags. If the ignore file is also in the todo list, then it is
   removed from the todo list. Creates nseq (100) inverse-folded seqs
   with a timeout of 10hrs using 5 cores. Note the inverse-fold uses 2
   cores/seq so :ncore 5 uses 10 cores total."

  [& args]
  (let [parse (fn [s] (-> (str/split #" " s) vec))
        [opts _ usage] (cli args
                            ["-t" "--todo" "files todo" :parse-fn parse
                             :default todo-files]
                            ["-d" "--done" "files done" :parse-fn parse
                             :default done-files]
                            ["-n" "--nseqs" "number of inverse seqs to create"
                             :parse-fn #(Integer/parseInt %) :default 100]
                            ["-nc" "--ncpu" "number of cpus to use"
                             :parse-fn #(Integer/parseInt %) :default 10]
                            ["-to" "--timeout" "timeout in hours"
                             :parse-fn #(Double/parseDouble %) :default 10]
                            ["-h" "--help" "usage" :default nil :flag true])
        {todo :todo done :done nseqs :nseqs} opts]
    (cond
     (opts :help)
     (print usage)

     args
     (doall (driver-create-inv todo done nseqs (opts :timeout)
                               :units :hr
                               :ncore (/ (opts :ncpu) 2)))

     :else
     (print usage))))

