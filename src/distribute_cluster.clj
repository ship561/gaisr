(ns distribute-cluster
  (:require [edu.bc.fs :as fs]
            [clojure.string :as str]
            [clojure.contrib.io :as io]
            [clojure.test :as test]
            [clojure-csv.core :as csv])
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.sequtils.files
        [edu.bc.bio.sequtils.snippets-files :only [read-clj]]
        edu.bc.utils))

(def ^{:private true} workdir (fs/join (fs/homedir) "bin/gaisr/trainset3/shuffled/"))

(defn- partition-jobs
  "partitions the files into groups so that all groups are similar in
  the number of sequences they contain"
  
  [file-list nbins]
  (let [file-lines (map (fn [f]
                          [(fs/basename f) (count (read-seqs f :info :names))]) 
                        file-list)]
    (prn file-lines)
    (loop [lim 30
           bins (reduce (fn [m [n files]]
                          (assoc m (* n (count files)) files))
                        {} (group-by second file-lines))]
      (cond
       (and (> (count bins) nbins)
            (pos? lim))
       (let [mk (apply min-key key bins) ;smallest key
             mk2 (apply min-key key (dissoc bins (first mk)));2nd smallest key
             new-bins (dissoc bins -1 (first mk) (first mk2))
             new-val (concat (second mk) (second mk2))
             new-key (apply + (map second new-val))]
         ;;if new-key exists then new-key=-1 and is dealt with
         ;;next round
         (recur (dec lim)
                (assoc new-bins (if (contains? new-bins new-key) -1 new-key) new-val)))

       #_(and (< (count bins) nbins)
            (pos? lim))
       #_(let [mk (apply max-key key bins) ;smallest key
             new-bins (dissoc bins Long/MAX_VALUE (first mk))
             [new-val1 new-val2] (partition-all (/ (count (second mk)) 3) (second mk))
             new-key1 (apply + (map second new-val1))
             new-key2 (apply + (map second new-val2))
             [new-key1 new-key2] (if (= new-key1 new-key2)
                                   [Integer/MAX_VALUE (dec Integer/MAX_VALUE)]
                                   [new-key1 new-key2])]
         ;;if new-key exists then new-key=-1 and is dealt with
         ;;next round
         (recur (dec lim)
                (assoc new-bins (if (contains? new-bins new-key1) Long/MAX_VALUE new-key1) new-val1
                       (if (contains? new-bins new-key2) (dec Long/MAX_VALUE) new-key2) new-val2)))
       :else
       bins))))


(defn- create-file-list
  "create the file list for putting into qsub file"
  
  [bins]
  (reduce (fn [v seqbin]
            (->> (map first seqbin)
                 (str/join " " )
                 (conj v)))
          [] (vals bins)))

(defn- command
  "Function generates a string to use as a command when running the
  function as a command line argument. any optional arguments such as
  \" -dfn distfn\" can by added."

  [function infiles workdir outfile & args]
  (apply str
         (str "lein run -m " function
              " -f " "\"" infiles "\""
              " -di " workdir
              " -o " outfile
              " -nc 16")
         (map #(str " " % " ") args)))

(defn- pbs-template
  "Currently used as a template for producing the pbs files for use on
  the shuffled data set in trainset3/shuffled"

  [outpbs cmdstr]
  (let [template {:shell "#!/bin/bash"
                  :resource "#PBS -l mem=5gb,nodes=1:ppn=16,walltime=100:00:00"
                  :working-dir "#PBS -d /home/peis/bin/gaisr/"
                  :command cmdstr}]
    (io/with-out-writer outpbs
      (println (template :shell))
      (println (template :resource))
      (println (template :working-dir) "\n")
      (println (template :command)))
    outpbs))

(defn combin-output [prefix out-file]
  (let [data (->> (str prefix ".\\d.out")
                  re-pattern 
                  (fs/re-directory-files "/home/peis/bin/gaisr/" )
                  (remove fs/empty? )
                  (mapcat read-clj ))]
    (io/with-out-writer out-file (prn (vec data)))))

(defn main-create-pbs [&  args]
  (let [[opts _ usage] (cli args
                            ["-f" "--file" "files to distribute evenly for jobs"
                             :parse-fn #(str/split % #" ") ;create list of files
                             :default nil]
                            ["-wdir" "--workdir" "dir in which files are located"
                             :default nil]
                            ["-di" "--pdir" "dir in which to produce pbs files"
                             :default nil]
                            ["-p" "--pbs" "prefix for pbs files to produce"
                             :default nil]
                            ["-o" "--out" "prefix for out files"
                             :default nil]
                            ["-i" "--start" "start number for pbs file suffix"
                             :parse-fn #(Integer/parseInt %)
                             :default 0]
                            ["-fn" "--function" "function to call in code"
                             :default "robustness/main-subopt-overlap"]
                            ["-efn" "--extrafn" "function used as arguments in command ie
                                                \"-efn subopt-overlap-neighbors\""
                             :default nil]
                            ["-n" "--partition-number" "number of partitions to make"
                             :parse-fn #(Integer/parseInt %)
                             :default 10])
        work-files (if (opts :workdir)
                     (map #(fs/join (opts :workdir) %) (opts :file))
                     (opts :file))
        parts (-> (partition-jobs work-files (opts :partition-number))
                  create-file-list)
        _ (prn :work parts)
        _ (cond
           (nil? args) (println usage)
           (nil? (opts :pbs)) (prn "need prefix for pbs files")
           (nil? (opts :out)) (prn "need prefix for out files")
           (nil? (opts :pdir)) (prn "need output directory for pbs files")
           (not (fs/exists? (opts :pdir))) (fs/mkdir (opts :pdir)))
        pbs-to-submit (map (fn [i j flist]
                             (pbs-template (fs/join (opts :pdir)
                                                    (str/join "." [(opts :pbs) i "pbs"]))
                                           (command (opts :function)
                                                    flist
                                                    (if (opts :workdir)
                                                      (opts :workdir)
                                                      (-> (opts :file) first fs/dirname))
                                                    (str (opts :out) "." j ".out")
                                                    (when (opts :extrafn) (opts :extrafn)))))
                           (iterate inc (opts :start))
                           (iterate inc 0)
                           parts)]
    pbs-to-submit))

(test/deftest proper-partitions
  (test/is (= (sum
               (map (fn [f]
                      (count (read-seqs f :info :names))) 
                    (fs/re-directory-files workdir "-NC*sto")))
              (->> (partition-jobs (fs/re-directory-files workdir "-NC*sto") 7)
                   keys
                   sum))))

(test/deftest key-val-same
  (test/is (= (->> (partition-jobs (fs/re-directory-files workdir "-NC*sto") 7)
                   keys)
              (->> (partition-jobs (fs/re-directory-files workdir "-NC*sto") 7)
                   vals
                   (map
                    #(reduce (fn [v [_ n]]
                               (+ v n))
                             0 %))))))
