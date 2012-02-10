(ns edu.bc.bio.zscore
  (:use edu.bc.utils
        net.n01se.clojure-jna)
  (:require [clojure.contrib.io :as io]
	    [edu.bc.fs :as fs]
            [clojure.contrib.string :as str]
            [incanter.stats :as stats]
            [clojure.java.shell :as shell])
  (:import (com.sun.jna Native)))

(defn jna-malloc
  "Create a 'C' level USB buffer of size SIZE.  Returns a pair (as a
   vector): [ptr-to-the-buffer the-buffer] Where ptr-to-the-buffer is
   a JNA/C ptr object and the-buffer is a java.nio.DirectByteBuffer
   object." [size]
  (let [buf (make-cbuf size)
        ptr (pointer buf)]
    [ptr buf]))

(def fdir "/home/peis/bin/zscore3/")

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
  (try
    (let [freq (frequencies s)
          gcat-per (double (/ (+ (get freq \g 0) (get freq \c 0)) length))
          at-per (double (/ (get freq \a 0) (+ (get freq \a 0) (get freq \u 0))))
          gc-per (double (/ (get freq \g 0) (+ (get freq \g 0) (get freq \c 0))))]
      (and (between? at at-per) (between? gc gc-per)
           (between? gcat gcat-per) (= length (count s))))
    (catch Exception e
      false)))

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
      (println ">" (str len "," gcat "," at "," gc))
      (println x))))

(defn generate-set-file [& {l :l :or {l (range 50 401 50)}}]
  (pmap (fn [length]
          (let [out-file (str fdir "set-" length ".txt")]
            (fs/mkdirs (fs/dirname out-file))
            (io/with-out-writer out-file (generate-set length 1000))))
        l))


(defn read-file2 [file]
  (doseq [[p s] (partition 2 (io/read-lines file))]
    (let [[l gcat at gc] (map #(Double/parseDouble %) (str/split #"," (subs p 2)))
          freq (frequencies seq)
          gcat-per (double (/ (+ (get freq \g 0) (get freq \c 0)) l))
          at-per (double (/ (get freq \a 0) (+ (get freq \a 0) (get freq \u 0))))
          gc-per (double (/ (get freq \g 0) (+ (get freq \g 0) (get freq \c 0))))]
      (when-not (and (between? at at-per) (between? gc gc-per)
		     (between? gcat gcat-per) (= l (count s)))
        (prn [p [l (count s)] [at at-per] [gc gc-per] [gcat gcat-per]])))))

(defn energy-of-seq [s]
  (jna-invoke Void RNA/read_parameter_file "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par")
  (jna-invoke Integer RNA/set_ribo_switch 1)
  (jna-invoke Void RNA/update_fold_params)
  (let [[ptr buf] (jna-malloc (inc (count s)))
        e (jna-invoke Float RNA/fold s ptr)
        ]
    e
    ;;(println "E =" e)
    ;; (println "sequence  =" i)
    ;;(println "Structure =" (.getString ptr 0 false))
    ))

(defn calc-energy-file [file]
  (reduce (fn [m x]
            (let [[p s] (apply map vector x)
                  es (map #(energy-of-seq %) s)]
              (assoc m (first p) [(stats/mean es) (stats/sd es)])))
          {} (group-by first (partition 2 (io/read-lines file)))))

(defn cef [file]
  (doseq [x (group-by first (partition 2 (io/read-lines file)))]
    (let [[p s] (apply map vector x)
          es (map #(energy-of-seq %) s)]
      (println (first p) [(stats/mean es) (stats/sd es)])))
   )

(defn calc-energy [& {l :l :or {l (range 50 401 50)}}]
  (pmap (fn [length]
          (let [in-file (str fdir "set-" length ".txt")
                out-file (str fdir "calc-set-" length ".txt")]
            (fs/mkdirs (fs/dirname out-file))
            (io/with-out-writer out-file (-> in-file
                                             ((fn [x] (shell/sh "cat" x)))
                                             ((fn [x] (shell/sh "RNAfold" "--noPS" "-P" "/usr/local/ViennaRNA-2.0.0/rna_andronescu2007.par" :in (get x :out))))
                                             ((fn [x] (print (get x :out))))))))
        l))

(defn map-energy 
  "takes a file and parses for the mfe to create a map where
   k=paramters and v=[(avg mfe of seq with k parameters) (sd mfe of seq
   with k parameters)]."
  [file]
  (reduce (fn [m [name mfe-list]]
            (let [es (map (fn [[_ _ st]]
                            (Double/parseDouble
                             (re-find #"\-*\d*.\d+" st)))
                          mfe-list)]
            (assoc m name [(stats/mean es) (stats/sd es) es])))
          {} (group-by first (partition 3 (io/read-lines file)))))

;; (io/with-out-writer (str fdir "sd.csv")
;;   (doseq [x (pmap (fn [in-file]
;;                     (map-energy in-file))
;;                   (for [i (range 50 401 50)]
;;                     (str fdir "calc-set-" i ".txt")))
;;           [att [_ sd _]] x]
;;     (println (subs att 2) "," sd)))

;; (do-text-file ["/home/peis/bin/overnight-set.txt"]
;;    :x)

;;(def runovernight
;;  (future (io/with-out-writer "/home/peis/bin/overnight-pmap.txt"
;;           (doseq [i (stat)] (doseq [j i] (prn j))))))

;; (def runovernight2
;;  (future (generate-set-file)))
