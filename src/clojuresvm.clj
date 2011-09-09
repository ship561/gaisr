(ns clojuresvm
  (:import (libsvm svm
                   svm_model
                   svm_problem
                   svm_node
                   svm_parameter)))

;;macro from
;;http://lilyx.net/2011/07/02/using-svm-support-vector-machine-from-clojure/
;;to set instance variables in java objects
(defmacro set-all! [obj m]
    `(do ~@(map (fn [e] `(set! (. ~obj ~(key e)) ~(val e))) m) ~obj))

(defn svm_node-construct [index value]
  (let [node (new svm_node)]
    (set-all! node {index index
                    value value})))

(defn svm_problem-construct [l y x]
  (let [prob (new svm_problem)]
    (set-all! prob {l l
                    y y
                    x x})))

(defn svm_parameter-construct [num_features options]
  (let [param (new svm_parameter)
        st {:C_SVC 0
            :NU_SVC 1
            :ONE_CLASS 2
            :EPSILON_SVR 3
            :NU_SVR 4}
        kt {:LINEAR 0
            :POLY 1
            :RBF 2
            :SIGMOID 3
            :PRECOMPUTED 4}
        m {:svm_type (get st (get options :svm_type) (st :C_SVC))
                 :kernel_type (get kt (get options :kernel_type) (kt :RBF))
                 :degree (get options :degree 3)
                 :gamma (get options :gamma (/ num_features))
                 :coef0 (get options :coef0 0)
                 :nu (get options :nu 0.5)
                 :cache_size (get options :cache_size 100)
                 :C (get options :C 1)
                 :eps (get options :eps 1e-3)
                 :p (get options :p 0.1)
                 :shrinking (get options :shrinking 1)
                 :probability (get options :probability 0)
                 :nr_weight (get options :nr_weight 0)
                 :weight_label (get options :weight_label (int-array 0))
                 :weight (get options :weight (double-array 0))}]
    (set-all! param {svm_type (m :svm_type)
                     kernel_type (m :kernel_type)
                     degree (m :degree)
                     gamma (m :gamma)
                     coef0 (m :coef0)
                     nu (m :nu)
                     cache_size (m :cache_size)
                     C (m :C)
                     eps (m :eps)
                     p (m :p)
                     shrinking (m :shrinking)
                     probability (m :probability)
                     nr_weight (m :nr_weight)
                     weight_label (m :weight_label)
                     weight (m :weight)})))

(defn addexample [row] ;;row=a vector of index value pairs
  (into-array svm_node (map (fn [[index value]]
                              (svm_node-construct index value))
                            row)))

(defn svm_node-toarray [m]
  (into-array (map (fn [row] 
                     (addexample row))
                   m)))

(defn svmlearn [labels train-set & {:as options}]
  (let [l (double-array labels)
        ts (svm_node-toarray train-set)
        param (svm_parameter-construct (count (first ts)) options)
        problem (svm_problem-construct (count l) l ts)]
    (if (nil? (svm/svm_check_parameter problem param))
      (svm/svm_train problem param)
      (prn "PROBLEM:" (svm/svm_check_parameter problem param)))))
      ;;(prn (svm/svm_predict model t))))))

(defn svmsavemodel [file model]
  (svm/svm_save_model file model))

(defn svmloadmodel [file]
  (svm/svm_load_model file))

(defn svmpredict [model test-set]
  (for [ts (svm_node-toarray test-set)]
    (do
      ;;(doseq [i ts] (println (.index i) (.value i) ""))
      (svm/svm_predict model ts))))
