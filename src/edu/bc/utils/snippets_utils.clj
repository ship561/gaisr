(ns edu.bc.utils.snippets-utils)

(defmacro with-out-appender [f & body]
  `(with-open [w# (clojure.java.io/writer ~f :append true)]
     (binding [*out* w#] ; I forgot this bit before.
       ~@body)))

(defn lazy-file-lines [file]
  (letfn [(helper [rdr]
            (lazy-seq
             (if-let [line (.readLine rdr)]
               (cons line (helper rdr))
               (do (.close rdr) nil))))]
    (helper (clojure.java.io/reader file))))

(defmacro print-proper
    "Do body in the repl normally. if a file f and body is provided, then
    it will do the body to the file."
    
    ([body] body)
    ([f body]
       `(with-open [w# (clojure.java.io/writer ~f)]
          (binding [*out* w#]
            ~body))))

(defmacro unless [pred a b]
  `(if (not ~pred) ~a ~b))
