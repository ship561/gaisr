(ns distribute-cluster
  (:require [edu.bc.fs :as fs]
            [clojure.string :as str]
            [clojure.contrib.io :as io]
            [clojure.test :as test])
  (:use [clojure.tools.cli :only [cli]]
        edu.bc.bio.sequtils.files
        edu.bc.utils))

(def ^{:private true} workdir (fs/join (fs/homedir) "bin/gaisr/trainset3/shuffled/"))

(defn- partition-jobs
  "partitions the files into groups so that all groups are similar in size"
  
  [workdir nbins]
  (let [file-lines (map (fn [f]
                          [(fs/basename f) (count (read-seqs f :info :names))]) 
                        (fs/re-directory-files workdir "-NC*sto"))]
    (loop [lim 30
           bins (group-by second file-lines)]
      (if (and (> (count bins) nbins)
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
        bins))))


(defn- create-file-list
  "create the file list for putting into qsub file"
  
  [bins]
  (reduce (fn [v seqbin]
            (->> (map first seqbin)
                 (str/join " " )
                 (conj v)))
          [] (vals bins)))

(defn- pbs-template
  "Currently used as a template for producing the pbs files for use on
  the shuffled data set in trainset3/shuffled"

  [outpbs infiles outfile]
  (let [template {:shell "#!/bin/bash"
                  :resource "#PBS -l mem=5gb,nodes=1:ppn=16,walltime=100:00:00"
                  :working-dir "#PBS -d /home/peis/bin/gaisr/"
                  :command (str "lein run -m robustness/main-subopt-overlap -f " infiles
                                " -di " workdir
                                " -o " outfile " -nc 16")}]
    (io/with-out-writer outpbs
      (println (template :shell))
      (println (template :resource))
      (println (template :working-dir) "\n")
      (println (template :command)))
    outpbs))

(defn main-create-pbs [&  args]
  (let [[opts _ usage] (cli args
                            ["-di" "--dir" "dir in which to produce pbs files"
                             :default nil]
                            ["-p" "--pbs" "prefix for pbs files to produce"
                                  :default nil]
                            ["-i" "--start" "start number for pbs file suffix"
                             :parse-fn #(Integer/parseInt %)
                             :default 0]
                            ["-n" "--partition-number" "number of partitions to make"
                             :parse-fn #(Integer/parseInt %)
                             :default 10])
        parts (-> (partition-jobs workdir (opts :partition-number))
                  create-file-list)
        _ (cond
           (nil? args) (println usage)
           (nil? (opts :pbs)) (prn "need prefix for pbs files")
           (nil? (opts :dir)) (prn "need output directory for pbs files")
           (not (fs/exists? (opts :dir))) (fs/mkdir (opts :dir)))
        pbs-to-submit (map (fn [i j flist]
                             (pbs-template (fs/join (opts :dir)
                                                    (str/join "." [(opts :pbs) i "pbs"]))
                                           flist
                                           (str "out." j ".out")))
                           (iterate inc (opts :start))
                           (iterate inc 0)
                           parts)]
    pbs-to-submit))

(test/deftest proper-partitions
  (test/is (= (sum
               (map (fn [f]
                      (count (read-seqs f :info :names))) 
                    (fs/re-directory-files workdir "-NC*sto")))
              (->> (partition-jobs workdir 7)
                   keys
                   sum))))

(test/deftest key-val-same
  (test/is (= (->> (partition-jobs workdir 7)
                   keys)
              (->> (partition-jobs workdir 7)
                   vals
                   (map
                    #(reduce (fn [v [_ n]]
                               (+ v n))
                             0 %))))))