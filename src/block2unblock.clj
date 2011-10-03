(ns stockholm
  (:require [clojure.contrib.io :as io])
  (:require [clojure.contrib.string :as str]))

(defn read-sto [f]
  (reduce (fn [m x] 
            (let [[name s] (str/split #"\s{2,}+" x)
                  s (str (get m name "") s)]
              (assoc m name s)))
          {} (io/read-lines f)))

(defn write-file [out-name in-file]
  (let [order (distinct (reduce (fn [m x] 
                                  (let [[name s] (str/split #"\s{2,}+" x)]
                                    (conj m name)))
                                [] (io/read-lines in-file)))
        sto (read-sto in-file)]
    (io/with-out-writer out-name 
      (doseq [i order]
        (if (< (count i) 15)
          (println i "\t\t" (sto i))
          (println i "\t" (sto i)))))))

