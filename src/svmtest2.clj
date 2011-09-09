(ns stuff
  (:require [clojure.contrib.io :as io]
            [clojuresvm :as csvm]
            [clojure.contrib.string :as str]))

(defn parse-line [line]
  (let [l (str/split #" " line)
        label (if (re-find #"\+" (first l))
                (subs (first l) 1)
                (first l))
        params (reduce (fn [m x]
                         (let [[index value] (str/split #":" x)]
                           (conj m [(Integer/parseInt index) (Double/parseDouble value)])))
                       [] (rest l))]
    [label params]))

(defn read-file [f]
  (partition 2 (reduce (fn [m ln]
                        (let [[label params] (parse-line ln)]
                          (conj m (Double/parseDouble label) params)))
                      [] (io/read-lines f))))

(defn train [file]
  (let [labels (map #(first %1) (read-file file))
        train-set (map #(second %1) (read-file file))]
    (csvm/svmlearn labels train-set :kernel_type :LINEAR))))

(defn predict [test-file model]
  (let [labels (map #(first %1) (read-file test-file))
        test-set (map #(second %1) (read-file test-file))]
    ;;(prn (count labels))
    (println "Accuracy:" (/ (count (filter true?
                                           (map #(= %1 %2) labels (csvm/svmpredict model test-set))))
                            (count labels)))
    (println "Accuracy to command line:" (/ (count (map #(= (Double/parseDouble %1) %2) (io/read-lines "/home/kitia/bin/libsvm-3.1/out.txt") (csvm/svmpredict model test-set)))
                                            (count (io/read-lines "/home/kitia/bin/libsvm-3.1/out.txt"))))))
