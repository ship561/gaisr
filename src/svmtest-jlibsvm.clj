(ns jlib
  (:import (java.util HashSet Vector))
  (:import (edu.berkeley.compbio.jlibsvm.kernel LinearKernel))
  (:import (edu.berkeley.compbio.jlibsvm ImmutableSvmParameterGrid))
  (:import (edu.berkeley.compbio.jlibsvm.binary C_SVC MutableBinaryClassificationProblemImpl))
  (:import (edu.berkeley.compbio.jlibsvm.util SparseVector)))
 
(defmacro set-all! [obj m]
    `(do ~@(map (fn [e] `(set! (. ~obj ~(key e)) ~(val e))) m) ~obj))
 
(defn into-sparsevec [m]
  (let [sv (new SparseVector (count m))
        sm (sort-by first m)]
    (set-all! sv {indexes (int-array (map first sm))
                  values (float-array (map second sm))})
    sv)
  )

(def svm (new C_SVC))
(def builder (ImmutableSvmParameterGrid/builder))
 
(set-all! builder {eps 1.0e-3
                   Cset (doto (new HashSet) (.add (float 1.0)))
                   kernelSet (doto (new HashSet) (.add (new LinearKernel)))})
 
(def param (. builder build))

(def x1 (into-sparsevec {1 1.0}))
(def x2 (into-sparsevec {1 -1.0}))
(def vx (new Vector [x1 x2]))
(def vy (new Vector [1 -1]))
 
(def prob (new MutableBinaryClassificationProblemImpl String (count vy)))
(doseq [x (map list vx vy)]
  (. prob (addExample (first x) (second x)))
  )

(def model (. svm (train prob param)))

(println (. model (predictLabel x1)))
(println (. model (predictLabel x2)))
