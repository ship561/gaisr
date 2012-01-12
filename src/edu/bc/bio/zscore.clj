(ns zscore
  (:use edu.bc.utils)
  (:require [clojure.contrib.io :as io]
	    [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]))

(def fdir "/home/peis/bin/zscore2/")

(defmacro <=_< [a b c]
  `(let [a# ~a
         b# ~b
         c# ~c]
     (and (<= a# b#) (< b# c#))))

(defn between? [theory actual]
  (< (- theory 0.025) actual (+ 0.025 theory)))

(defn next-base [gcat at gc]
  (let [r (rand)]
    ;;(pr r (* gc gcat) (+ gcat (* at (- 1 gcat))))
    (cond
     (< r (* gc gcat)) "g"
     (<=_< (* gc gcat) r gcat) "c"
     (<=_< gcat r (+ gcat (* at (- 1 gcat)))) "a"
     :else "u")))

(defn valid-seq?
  [length gcat at gc s]
  (let [freq (frequencies s)
        gcat-per (double (/ (+ (get freq \g 0) (get freq \c 0)) length))
        at-per (double (/ (get freq \a 0) (+ (get freq \a 0) (get freq \u 0))))
        gc-per (double (/ (get freq \g 0) (+ (get freq \g 0) (get freq \c 0))))]
    (and (between? at at-per) (between? gc gc-per)
         (between? gcat gcat-per) (= length (count s)))))

(defn generate-seqs [length gcat at gc]
  (let [generate-seq (fn [] (seq (repeatedly length #(next-base gcat at gc))))]
    (distinct (seque (filter (fn [s]
                               (valid-seq? length gcat at gc s))
                             (repeatedly #(str/join "" (generate-seq))))))))

(defn generate-set [len nseq]
  (doseq [gcat (range 0.25 0.76 0.05)
          at (range 0.25 0.76 0.05)
          gc (range 0.25 0.76 0.05)]
    (doseq [x (take nseq (generate-seqs len gcat at gc))]
      (println ">" [len gcat at gc])
      (println x))))

;; (defn generate-set-files [dir]
;;   (pmap #(generate-set-file %1 (str dir "seq-data-" % ".txt"))
;;         (range 50 401 50)))

(defn generate-set-file []
  (pmap (fn [length]
          (let [out-file (str fdir "set-" length ".txt")]
            (fs/mkdirs (fs/dirname out-file))
            (io/with-out-writer out-file (generate-set-file length 1000))))
        (range 50 401 50)))


(defn read-file2 [file]
  (doseq [[p s] (partition 2 (io/read-lines file))]
    (let [[l gcat at gc] (map #(Double/parseDouble %) (str/split #"," (subs p 2)))
          seq (drop 1 (str/split #"" s))
          freq (frequencies seq)
          gcat-per (double (/ (+ (get freq "g" 0) (get freq "c" 0)) l))
          at-per (double (/ (get freq "a" 0) (+ (get freq "a" 0) (get freq "u" 0))))
          gc-per (double (/ (get freq "g" 0) (+ (get freq "g" 0) (get freq "c" 0))))]
      (when-not (and (between? at at-per) (between? gc gc-per)
		     (between? gcat gcat-per) (= l (count s)))
        (prn [p [l (count s)] [at at-per] [gc gc-per] [gcat gcat-per]])))))




;; (do-text-file ["/home/peis/bin/overnight-set.txt"]
;;    :x)



;; (defn ratio [m k1 k2]
;;   (double (/ (get m k1 0) (+ (get m k1 0) (get m k2 0)))))

;; (defn stat [in-set]
;;   (map (fn [x]
;;          (map (fn [[param y]]
;;                 (let [fr (frequencies (flatten y))
;;                       total (reduce (fn [s [k v]]
;;                                       (+ s v))
;;                                     0 fr)
;;                       fr-per (reduce (fn [m [k v]]
;;                                        (assoc m k (double (/ v total))))
;;                                      {} fr)
;;                       gcat (double (+ (get fr-per "g" 0) (get fr-per "c" 0)))
;;                       at (ratio fr-per "a" "t")
;;                       gc (ratio fr-per "g" "c")
;;                       fr-per (assoc fr-per "gcat" gcat "at" at "gc" gc)]
;;                   [param (sort-by first fr-per) total]))
;;               x))
;;        in-set))

;; (defn print-stat [in-set]
;;   (doseq [i (stat in-set)]
;;     (doseq [j i]
;;       (let [[p per n] j
;;             gcat (nth p 1)
;;             at (nth p 2)
;;             gc (nth p 3)
;;             at-per (second (nth per 1))
;;             gc-per (second (nth per 4))
;;             gcat-per (second (nth per 5))
;;             between? (fn [theory actual] (< (* 0.99 theory) actual (* 1.01 theory))) ]
;;         (when-not (or (between? at at-per) (between? gc gc-per) (between? gcat gcat-per))
;;           (prn p ", at" [at at-per] ", gc" [gc gc-per] ", gcat" [gcat gcat-per]))))))


;;(def runovernight
;;  (future (io/with-out-writer "/home/peis/bin/overnight-pmap.txt"
;;           (doseq [i (stat)] (doseq [j i] (prn j))))))

;; (def runovernight2
;;  (future (io/with-out-writer "/home/peis/bin/overnight-set.txt"
;;            (print-set (generate-set)))))
