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
                        file-list)
        bins (map (fn [x] [(second x) [x]]) (take nbins file-lines))]
    (loop [current-bins (sort-by first bins)
           remaining-files (drop nbins file-lines)]
      (if (seq remaining-files)
        (let [i (first remaining-files)
              j (second (first current-bins))]
          ;;(prn :remaining i :cbins j) 
          (recur (let [update-bin (conj j i)
                       cnt (sum (map second update-bin))]
                   ;;(prn :update update-bin :cnt cnt)
                   (->> [cnt update-bin]
                        (conj (rest current-bins))
                        (sort-by first )))
                 (rest remaining-files)))
        current-bins))))


(defn- create-file-list
  "create the file list for putting into qsub file"
  
  [bins]
  (reduce (fn [v seqbin]
            (->> (second seqbin)
                 (map first)
                 (str/join " " )
                 (conj v)))
          [] bins))

(defn- command
  "Function generates a string to use as a command when running the
  function as a command line argument. any optional arguments such as
  \" -dfn distfn\" can by added."

  [function infiles workdir outfile ncore & args]
  (apply str
         (str "lein run -m " function
              " -f " "\"" infiles "\""
              " -di " workdir
              " -o " outfile
              " -nc " ncore)
         (map #(str " " % " ") args)))

(defn- pbs-template
  "Currently used as a template for producing the pbs files for use on
  the shuffled data set in trainset3/shuffled"

  [outpbs cmdstr ncore]
  (let [template {:shell "#!/bin/bash"
                  :resource (str "#PBS -l mem=5gb,nodes=1:ppn=" ncore ",walltime=100:00:00")
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
                  (fs/re-directory-files (fs/join (fs/homedir) "bin/gaisr"))
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
                             :default "robustness.neutrality/main-subopt-overlap"]
                            ["-efn" "--extrafn" "function used as arguments in command ie \"-dfn subopt-overlap-neighbors\""
                             :default nil]
                            ["-n" "--partition-number" "number of partitions to make"
                             :parse-fn #(Integer/parseInt %)
                             :default 10]
                            ["-nc" "--ncore" "number of cores to use per partition"
                             :parse-fn #(Integer/parseInt %)
                             :default 16])
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
                                                    (opts :ncore)
                                                    (when (opts :extrafn) (opts :extrafn)))
                                           (opts :ncore)))
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
