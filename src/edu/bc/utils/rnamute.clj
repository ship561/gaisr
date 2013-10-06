(ns edu.bc.utils.rnamute
  (:require [clojure.string :as str]
            [clojure.contrib.io :as io]
            [clojure.java.shell :as shell]
            [clojure.test :as test]
            [edu.bc.fs :as fs]
            )
  (:use [clojure.tools.cli :only [cli]]
        [edu.bc.bio.sequtils.files :only (read-seqs)]
        [robustness :only [valid-seq-struct]]
        [snippets-analysis :only [alignment-quality]])
  )

(def mute-dir (str (fs/homedir) "/bin/RNAMute/"))

(defn seq->file
  "Takes a sequence and prints it to a temp file in the RNAmute directory."

  [inseq]
  (let [tmp (fs/tempfile "-seq-" "" mute-dir)]
    (io/with-out-writer tmp (println inseq))
    tmp))

(defn sto->seqs
  "Takes a sto and prints each of the sequences (gaps removed) to temp
  files in the RNAmute directory"

  [insto]
  (->> (map #(when (nil? (re-find #"[^ACGTUacgtu\.]" %))
               (-> (str/replace % #"\." "") seq->file))
            (read-seqs insto))
       (remove nil? )))

(defn run-rnamute
  "Shell call to RNAmute"

  [infile]
  (shell/sh "java" "-jar" "RNAmute.jar" infile :dir mute-dir))

(defn parse-output
  "Parses the output from RNAMute in the result_table.html file.
  Returns a vector [dot-bracket-distance shapiro-notation-distance].

  example file \"/home/kitia/Desktop/RNAMute/RESULT_TABLE.html\""

  []
  (let [out (-> (str mute-dir "RESULT_TABLE.html") io/read-lines second)
        parse-for (fn [re] (-> (re-find re out) second Double/parseDouble))
        dotdist (parse-for #"Dot-Bracket Distances is: (\-?\d*\.?\d+)")
        shapiro (parse-for #"Shapiro's Distances is: (\-?\d*\.?\d+)")]
    [dotdist shapiro]))

(defn RNAmute-dist
  "main function to wrap the RNAmute shell calls with the proper input
  files. Returns a list of vectors in the format [dot-bracket-distance
  shapiro-notation-distance]"

  [insto]
  (let [norm-dist (fn [d L] (- 1 (/ d L)))] ;1-d/L
    (for [s (sto->seqs insto)
          :let [len (-> (io/read-lines s)
                        doall
                        first
                        count)]]
      ;;first parses the RNAmute output then normalizes it in the mapfn
      (mapv #(norm-dist % len)
           (do (run-rnamute s)
               (parse-output))))))

(defn -main [& args]
  (let [[opts _ usage] (cli args
                            ["-f" "--file" "file(s) to check neutrality for"
                             :parse-fn #(str/split % #" ") ;create list of files
                             :default nil]
                            ["-o" "--outfile" "file to write to"
                             :default (str (fs/homedir) "/bin/gaisr/robustness/rnamute-data.clj")]
                            ["-di" "--dir" "dir in which files are located"
                             :default nil]
                            ["-d" "--debug" "debug program prints to repl"
                             :default nil
                             :flag true]
                            )
        files (if (opts :dir)
                (map #(fs/join (opts :dir) %) (opts :file))
                (opts :file))
        results (mapv #(do (prn %)
                           (conj (RNAmute-dist %) %))
                      files)
        outtxt (fn [] (doseq [[nm & vs] results
                             [bpdist shapirodist] vs]
                       (print (str nm "," bpdist ",bp-distance\n"
                                   nm "," shapirodist ",shapiro-dist\n"))))]
    (if (opts :debug)
      (do (prn results) (outtxt))
      (io/with-out-writer (opts :outfile)
        (println "#data for rnamute program running it on all seqs in trainset2/pos with re=\".\\d.sto$\". returns a vector of vectors of the average distance already normalized using <1-d/L>. Vector format is [dotbracket-dist shapiro-dist].")
        (outtxt)))))


(test/deftest- footest
  
  (test/is (= (RNAmute-dist (str (fs/homedir) "/bin/gaisr/trainset2/pos/RF00167-seed.3.sto"))
              '([0.827186274509804 0.7580490196078431] [0.8586 0.72377] [0.9147339999999999 0.84964])))
  (test/is (= (map RNAmute-dist (take 3 (fs/re-directory-files (str (fs/homedir) "/bin/gaisr/trainset2/pos/") #"RF00167.*\.\d\.sto")))
            '(([0.9276 0.88634] [0.8288529411764706 0.7038333333333333] [0.9322058823529412 0.8767549019607843] [0.8879313725490197 0.8069411764705883] [0.7681078431372549 0.6504313725490196] [0.85107 0.78687]) ([0.8617254901960785 0.7863725490196078] [0.8586 0.72377] [0.956067 0.9394] [0.879343137254902 0.8015196078431372] [0.6476530612244897 0.4582142857142857] [0.71634 0.64517]) ([0.9370117647058823 0.9017058823529411] [0.956067 0.9394] [0.8508921568627451 0.7728137254901961] [0.9096 0.8439] [0.9328470588235294 0.9161225490196079] [0.8569166666666667 0.7571979166666667])))))
