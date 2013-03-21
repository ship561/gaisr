(ns snippets-partition-function
  (:use clojure.contrib.map-utils))

(def ztable (atom {}))

(defn- bp? [b1 b2] ;;b1=base1 b2=base2
  (let [bp #{"AU" "UA" "GC"  "CG" "GU" "UG"}]
    (contains? bp (str b1 b2))))

(defn- e
  "Energy function = 1 to count structures
  but needs to be changed to actual energies."
  
  [i j S count]
  (let [s1 (subs S i (inc i))
        s2 (subs S j (inc j))
        R 2
        T 310
        E (fn [b1 b2] (if (bp? b1 b2) 1 0))] ;;Energy of basepair, E(basepair)
    (if count
      (E s1 s2)
      (Math/exp (/ (E s1 s2) R T -1)))))

(defn Z
  "Calculates the partition function for a structure (from sequence S)
  from i (starts at 0) to j. Functions based on Nussinov Z(1,n) =
  Z(1,n-1) + Zb(1,n) + sum(Z(1,k-1)*Zb(k,n) for 1<k<n) where Zb(1,n) =
  exp(-a(1,n)/RT)*Z(2,n-1) and Zb(k,n)=exp(-a(k,n)/RT)*Z(k+1,n-1)

  The probability of a base pair i,j occuring is given by:
  Pr(i,j) = (* Z(1,i-1) E(s1 s2) Z(i+1,j-1) Z(j+1,n))/Z(1,n)

  count= only count structures, u=min loop size"
  
  [i j S & {:keys [count u]
            :or {count false ;count structures
                 u 3}}] ;min loop size
  (let [_ (reset! ztable {})
        z (fn z [i j S]
            (let [S (.toUpperCase S)
                  n (- j i -1) 
                  result (if (<= (- j i) u) 1
                             (+' (lazy-get @ztable [i (dec j)] (z i (dec j) S)) ;;j unpaired
                                (*' (e i j S count)
                                   (lazy-get @ztable [(inc i) (dec j)] (z (inc i) (dec j) S))) ;;i,j pair
                                (reduce (fn [x k] ;;k,j paired for an intermediate a<k<b
                                          (+' x (*' (e k j S count)
                                                  (lazy-get @ztable [i (dec k)] (z i (dec k) S))
                                                  (lazy-get @ztable [(inc k) (dec j)] (z (inc k) (dec j) S)))))
                                        0 (range (inc i) (- j u)))))]
              (swap! ztable #(assoc % [i j] result))
              result))]
    (z i j S)
    (get @ztable [i j])))

(defn bp-prob
  "Pr(i,j) = (* Z(1,i-1) E(s1 s2) Z(i+1,j-1) Z(j+1, n))/(Z(1,n)"
  
  [i j s]
  (let [u 3
        n (-> s count dec)]
    (/ (* (Z 0 (dec i) s :count true :u u) ;Z(1,i-1)
          (if (<= (- j i) u) 0 (e i j s true)) ;energy of i,j bp
          (Z (inc i) (dec j) s :count true :u u) ;Z(i+1,j-1)
          (Z (inc j) n s :count true :u u)) ;Z(j+1,n)
       (Z 0 n s :count true :u u)))) ;partition function


