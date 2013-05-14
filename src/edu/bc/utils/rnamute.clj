(ns edu.bc.utils.rnamute
  (:require [clojure.string :as str]
            [clojure.contrib.io :as io]
            [clojure.java.shell :as shell]
            [edu.bc.fs :as fs])
  (:use [edu.bc.bio.sequtils.files :only (read-seqs)]
        ))

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
  (map #(-> (str/replace % #"\." "") seq->file) (read-seqs insto)))

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
        dotdist (parse-for #"Dot-Bracket Distances is: (\d+.?\d+)")
        shapiro (parse-for #"Shapiro's Distances is: (\d+.?\d+)")]
    [dotdist shapiro]))

(defn RNAmute-dist
  "main function to wrap the RNAmute shell calls with the proper input
  files. Returns a list of vectors in the format [dot-bracket-distance
  shapiro-notation-distance]"

  [insto]
  (let [norm-dist (fn [d L] (- 1 (/ d L)))]
    (for [s (sto->seqs insto)
          :let [len (-> (io/read-lines s)
                        first
                        count)]]
      (mapv #(norm-dist % len)
           (do (run-rnamute s)
               (parse-output))))))
