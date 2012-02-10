(ns libsvm2weka
  (:require [clojure.contrib.string :as str]
            [clojure.contrib.io :as io])
  (:import weka.classifiers.functions.LibSVM
           weka.classifiers.lazy.IBk
           weka.classifiers.Evaluation))

(defn txt2csv [in outfile]
  (io/with-out-writer outfile
    (println "length,gcat,at,gc,class")
    (doseq [i in]
      (let [fun (fn [x] (second (str/split #":" x)))
            [class a b c d] (str/split #" " i)]
        (println (str/join "," [(fun a) (fun b) (fun c) (fun d) class]))))))

;;processes the output from the classifiers and saves only the
;;predicted values
(defn getpredicted [output]
  (reduce (fn [p ln]
              (let [[instnum actual predicted error] (remove empty? (str/split #" " ln))]
                (conj p (Double/parseDouble predicted))))
          [] (drop 5 (str/split-lines output))))

;;works to predict based on test data set
(defn wekalibsvm [trainfile testfile modelname]
  (let [svm (new LibSVM)
        ;;trainfile "/home/kitia/mean1.csv"
        ;;testfile "/home/kitia/test.csv"
        ;; -K 2 -D 3 -G 0.0 -R 0.0 -N 0.5 -M 40.0 -E 0.0010 -P 0.1
        model (if (= "mean" modelname) "mean.model" "sd.model")
        options (into-array java.lang.String [;;"-t" trainfile 
                                              "-T" testfile
                                              ;;if using a model
                                              ;;can't use svm
                                              ;;options -Z or -S
                                              ;;"-d" "modeltosave.model"
                                              ;;"-l" "modeltoload.model"
                                              ;;"-G" "0.0"
                                              ;;"-C" "1.0"
                                              "-x" "5"
                                              "-p" "0"
                                              "-i"
                                              ;;"-S" "3" 
                                              ;;"-Z"
                                              "-l" model
                                              ])]
    (Evaluation/evaluateModel svm options)))

(defn wekaNN [trainfile testfile]
  (let [nn (new IBk)
        ;;trainfile "/home/kitia/mean1.csv"
        ;;testfile "/home/kitia/test.csv"
        options (into-array java.lang.String ["-t" trainfile 
                                              "-T" testfile
                                              ;;if using a model
                                              ;;can't use svm
                                              ;;options -Z or -S
                                              ;;"-d" "modeltosave.model"
                                              ;;"-l" "modeltoload.model"
                                              "-x" "5"
                                              "-p" "0"
                                              "-i"
                                              "-K" "10"])]
    (Evaluation/evaluateModel nn options)))
