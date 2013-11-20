(ns smith-waterman)

(defn- array-keys
  "positions of the array to work on"

  [s1 s2]
  (for [i (range (count s1)) ;initialize scoring array. similar to a sparse matrix
        j (range (count s2))]
    [i j]))

(defn- init-array
  "initial conditions of the array to fill in for the gapped row/col"
  
  [s1 s2 locations gap global]
  (reduce (fn [m [i j]]
            (assoc m [i j]
                   (cond
                    (= i j 0) [0 [0 0] "-" "-"]
                    (= i 0) [(if global (* gap j) 0) [0 (dec j)] "-" (.charAt s2 j)]
                    (= j 0) [(if global (* gap i) 0) [(dec i) 0] (.charAt s1 i) "-"])))
          {} (filter #(or (zero? (first %))
                          (zero? (second %))) locations)))

(defn- fill-array
  
  [locations s1 s2 match mis gap global]
  (reduce (fn [m [i j]];;score array format
            (let [maxa (fn [coll]
                         (->> (filter #(= (apply max (map second coll))
                                          (second %)) coll)
                              (sort-by first)
                              first))
                  d (first (get m [(dec i) (dec j)])) ;;score match/mismatch (diagonal)
                  u (first (get m [(dec i) j])) ;;score deletion (above)
                  l (first (get m [i (dec j)])) ;;score insertion (left)
                  aa1 (.charAt s1 i) ;;current char in s1
                  aa2 (.charAt s2 j) ;;current char in s2
                  [from score] ;;chooses from d, u, l and scores associated with it.
                  (maxa [(if (= aa1 aa2 ) 
                           [:d (+ d match)]
                           [:d (+ d mis)])
                         [:u (+ u gap)]
                         [:l (+ l gap)]])]
              (assoc m [i j] ;;insertion of the best score into the matrix
                     (case from
                       :d [score [(dec i) (dec j)] aa1 aa2]
                       :u [score [(dec i) j] aa1 "-"]
                       :l [score [i (dec j)] "-" aa2]))))
          (init-array s1 s2 locations gap global)
          (remove #(or (zero? (first %))
                       (zero? (second %))) locations)))

(defn- trace-back

  [score start-loc H]
  (loop [loc start-loc
         aln_s1 ""
         aln_s2 ""]
    (let [[_ [i j] a1 a2] (get H loc)] ;;stores the next location [score[i j] to go to in H]
      (if (not= [0 0] [i j])
        (recur [i j]
               (str a1 aln_s1) ;;builds strings up from the right to left
               (str a2 aln_s2))
        (if (= "-" a1 a2)
          [score aln_s1 aln_s2]
          [score (str a1 aln_s1) (str a2 aln_s2)])))))

(defn sw [seq1 seq2 & {:keys [global
                              match-weight
                              mismatch-weight
                              gap-weight]
                       :or {global false
                            match-weight 2
                            mismatch-weight -1
                            gap-weight -1}} ]
  (let [s1 (str "-" seq1)
        s2 (str "-" seq2)
        match match-weight  ;match
        mis mismatch-weight ;mismatch
        gap gap-weight ;gap 
        H (fill-array (array-keys s1 s2) s1 s2 match mis gap global) ;creates score matrix
        start (if global
                (get H (mapv dec [(count s1) (count s2)]))
                (-> (sort-by #(-> % second first) H) ;finds highest value in matrix
                    last
                    second))
        start-locs (map first (filter #(= start (val %)) H))] ;;starts traceback from this highest value
    (map #(trace-back (first start) % H) start-locs)))
