;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                    S E Q U T I L S . D I S T S                           ;;
;;                                                                          ;;
;;                                                                          ;;
;; Copyright (c) 2011-2012 Trustees of Boston College                       ;;
;;                                                                          ;;
;; Permission is hereby granted, free of charge, to any person obtaining    ;;
;; a copy of this software and associated documentation files (the          ;;
;; "Software"), to deal in the Software without restriction, including      ;;
;; without limitation the rights to use, copy, modify, merge, publish,      ;;
;; distribute, sublicense, and/or sell copies of the Software, and to       ;;
;; permit persons to whom the Software is furnished to do so, subject to    ;;
;; the following conditions:                                                ;;
;;                                                                          ;;
;; The above copyright notice and this permission notice shall be           ;;
;; included in all copies or substantial portions of the Software.          ;;
;;                                                                          ;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          ;;
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       ;;
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    ;;
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   ;;
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ;;
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    ;;
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          ;;
;;                                                                          ;;
;; Author: Jon Anthony                                                      ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns edu.bc.bio.sequtils.dists

  "Various sequence, sequence families, alignments and alignment
   families distributions.  Intended for use as \"background\"
   distributions in relative entropy and other such comparaitive or
   measure calculations.  These distributions are focused on microbial
   NC genomes and RFAM alignment families.  At least to start.
   Others, such as metagenomic results and NZ genomes could be here in
   future."

  (:require [clojure.contrib.string :as str]
            [clojure.contrib.str-utils :as stru]
            [clojure.set :as set]
            [clojure.contrib.seq :as seq]
            [clojure.zip :as zip]
            [clojure.contrib.io :as io]
            [clj-shell.shell :as sh]
            [edu.bc.fs :as fs])
  (:use clojure.contrib.math
        edu.bc.utils
        edu.bc.utils.probs-stats
        edu.bc.bio.seq-utils
        edu.bc.bio.sequtils.files
        edu.bc.bio.sequtils.info-theory
        [clojure.contrib.condition
         :only (raise handler-case *condition* print-stack-trace)]
        [clojure.contrib.pprint
         :only (cl-format compile-format)]
        ))



;;; ----------------------------------------------------------------------
;;;
;;; IUPAC (International Union of Pure and Applied Chemistry) codes
;;; for nucleotides and various "groups" of nucleotides.
;;;

(def +IUPAC+
     ^{:doc "These are the codes for \"bases\" used in alignments"}
     {\A "Adenine"
      \C "Cytosine"
      \G "Guanine"
      \T "Thymine"
      \U "Uracil"
      \R "AG"
      \Y "CUT"
      \S "GC"
      \W "AUT"
      \K "GUT"
      \M "AC"
      \B "CGUT"
      \D "AGUT"
      \H "ACUT"
      \V "ACG"
      \N "any"
      \. "gap"
      \- "gap"})

(def +NONSTD-RNA+
     ^{:doc "These are the codes for all non standard \"bases\" used
             in alignments"}

     {\T "Thymine"
      \R "AG"       ; purines
      \Y "CUT"      ; pyrimidines
      \S "GC"
      \W "AUT"
      \K "GUT"
      \M "AC"
      \B "CGUT"
      \D "AGUT"
      \H "ACUT"
      \V "ACG"
      \N "any"
      \. "gap"
      \- "gap"})


;;; ----------------------------------------------------------------------
;;;
;;; Some bacground frequency and probability distributions for all NC
;;; genomes and RFAM families.  May well be better to place in
;;; database and be pulled in and initialized at startup time...
;;;

(defn print-dist
  "Prints a distribution DIST, a frequency or probability distribution
   map, in sorted order of most occurances / highest probability to
   lowest.  One item per line.
  "
  [dist]
  (doseq [[k v] (sort-by second > dist)]
    (print " ") (pr k) (print " ") (pr v) (println \,)))


(defn prn-dist
  "Print readably the distribution DIST to file FILE.  This means that
   FILE can be used as a secondary cache and so can be simply read
   with Clojures reader to get the internal representation back."
  [dist file]
  (io/with-out-writer (fs/fullpath file)
    (prn dist)))



;;; (bg-freqs-probs 2 "/data2/BioData/GenomeSeqs"
;;;                 :ftypes [".fna"] :norm false :par 4)

(def +bg-nc-genomes+
     ^{:doc
       "Dinucleotide Frequencies for all NC microbial genomes.  In
        addition to the 16 \"standard\" 4 base pair (A,T,G,C)
        permutations, includes many non standard pairs (see
        +NONSTD-RNA+"}

     {"GC" 395069986,
      "CG" 377952477,
      "TT" 326305148,
      "AA" 326263282,
      "CC" 306796217,
      "GG" 306596839,
      "AT" 299033486,
      "CA" 291400632,
      "TG" 291273508,
      "TC" 287515297,
      "GA" 287493041,
      "AG" 251060051,
      "CT" 251035265,
      "AC" 237803474,
      "GT" 237723100,
      "TA" 209003358,
      "NN" 11849, "YT" 693, "NA" 688, "AN" 676, "TN" 664, "NT" 626, "NG" 480,
      "NC" 470, "GN" 468, "CN" 458, "AR" 403, "TY" 359, "CY" 335, "GY" 313,
      "RA" 295, "RG" 245, "RC" 229, "YG" 227, "AY" 227, "CR" 222, "CS" 189,
      "MC" 184, "CM" 177, "YC" 177, "SC" 172, "AK" 171, "RT" 170, "WA" 162,
      "TW" 156, "AM" 155, "GS" 155, "GR" 153, "TR" 150, "KG" 149, "YA" 143,
      "SG" 143, "MA" 140, "AS" 139, "GK" 135, "KA" 133, "TK" 127, "KC" 126,
      "ST" 125, "KT" 121, "MT" 110, "CK" 108, "WT" 107, "SA" 106, "GM" 102,
      "AW" 101, "CW" 99, "MG" 99, "TM" 93, "WG" 84, "TS" 76, "WC" 74, "GW" 68,
      "RR" 25, "MM" 18, "SS" 17, "KK" 14, "YY" 14, "YW" 14, "KR" 13, "SM" 12,
      "WW" 11, "SY" 11, "MY" 11, "WR" 11, "RK" 10, "WS" 10, "WY" 9, "RY" 9,
      "KY" 9, "YK" 9, "YR" 9, "SR" 9, "KM" 8, "SW" 8, "KS" 8, "MW" 8, "RM" 8,
      "BG" 7, "CH" 7, "MR" 7, "RW" 7, "GB" 7, "RS" 6, "VA" 6, "BT" 6, "CV" 6,
      "WM" 6, "SK" 6, "YS" 6, "CB" 6, "AD" 5, "MS" 5, "KW" 5, "VC" 5, "AV" 5,
      "HC" 5, "CD" 4, "TV" 4, "DG" 4, "HT" 4, "TB" 4, "GV" 4, "VG" 4, "YM" 4,
      "BA" 4, "DC" 4, "HH" 3, "GH" 3, "WK" 3, "SN" 3, "VT" 3, "AB" 2, "DT" 2,
      "TD" 2, "YN" 2, "DA" 2, "GD" 2, "NK" 2, "MK" 2, "HG" 2, "NM" 2, "BC" 1,
      "MN" 1, "VW" 1, "KN" 1, "VY" 1, "NR" 1, "AH" 1, "DM" 1, "MV" 1, "NW" 1,
      "HR" 1, "NY" 1, "YH" 1, "BR" 1, "TH" 1, "WN" 1, "HA" 1, "RN" 1, "VS" 1,
      "WV" 1})

(def +bg-nc-genomes-probs+
     ^{:doc
       "Dinucleotide probabilities for all NC microbial genomes."}
     (probs +bg-nc-genomes+))


(def +bg-nc-genomes-sym+
     ^{:doc
       "Dinucleotide Frequencies for all NC microbial genomes where
        pairs are treated symmetrically, i.e., XY is considerted = YX.
        In addition to the 10 \"standard\" 4 base pair (A,T,G,C)
        symmetric permutations, includes many non standard pairs (see
        +NONSTD-RNA+"}

     {"GC" 773022463,
      "GA" 538553092,
      "TC" 538550562,
      "CA" 529204106,
      "TG" 528996608,
      "AT" 508036844,
      "TT" 326305148,
      "AA" 326263282,
      "CC" 306796217,
      "GG" 306596839,
      "NN" 11849, "NA" 1364, "TN" 1290, "YT" 1052, "NG" 948, "NC" 928,
      "AR" 698, "GY" 540, "CY" 512, "RC" 451, "RG" 398, "AY" 370, "CS" 361,
      "MC" 361, "RT" 320, "AK" 304, "GS" 298, "AM" 295, "KG" 284, "TW" 263,
      "WA" 263, "TK" 248, "AS" 245, "KC" 234, "MT" 203, "ST" 201, "GM" 201,
      "CW" 173, "WG" 152, "RR" 25, "KR" 23, "YW" 23, "MM" 18, "RY" 18,
      "KY" 18, "WR" 18, "WS" 18, "SS" 17, "SY" 17, "SM" 17, "MY" 15, "RM" 15,
      "SR" 15, "KK" 14, "YY" 14, "BG" 14, "KS" 14, "MW" 14, "CH" 12, "WW" 11,
      "VA" 11, "CV" 11, "KM" 10, "BT" 10, "CD" 8, "KW" 8, "GV" 8, "TV" 7,
      "AD" 7, "CB" 7, "DG" 6, "BA" 6, "GH" 5, "HT" 5, "DT" 4, "HH" 3,
      "YN" 3, "SN" 3, "NK" 3, "NM" 3, "VW" 2, "NR" 2, "AH" 2, "NW" 2,
      "VY" 1, "DM" 1, "MV" 1, "HR" 1, "YH" 1, "BR" 1,
      "VS" 1})

(def +bg-nc-genomes-sym-probs+
     ^{:doc
       "Dinucleotide probabilities for symmetric frequencies for all
       NC microbial genomes."}
     (probs +bg-nc-genomes-sym+))


;;; (bg-freqs-probs 1 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam+
     (let [fs { "." 634856, "A" 157097, "G" 154623, "C" 130024, "U" 128626,
                "N" 7, "Y" 1, "S" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM" :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam-di+
     (let [fs { ".." 570773, "AA" 46303, "GG" 40160, "GA" 36784, "AG" 34170,
                "GC" 32701, "UU" 30849, "CC" 30216, "UG" 29804, "CG" 29361,
                "GU" 28772, "CA" 28589, "AU" 28537, "AC" 27831, "UC" 27120,
                "UA" 25212, "CU" 25150, ".G" 20125, "A." 18870, ".A" 18021,
                "C." 15578, "G." 15092, "U." 14374, ".U" 14005, ".C" 11403,
                "CN" 3, "N." 2, "NU" 2, "NA" 2, "UN" 2, "SU" 1, "UY" 1,
                "GN" 1, "AN" 1, "NC" 1, "YU" 1, "US" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM"
;;;                 :fsps-fn cc-combins-freqs-probs
;;;                 :cols true :ftypes [".sto"] :norm true :par 2)
(def +bg-freqs-probs-all-rfam-cols+
     (let [fs {".." 165988848, "GG" 15598223, "AA" 14862080, "CC" 11859618,
               "UU" 8978659, "AG" 6619716, "UA" 6160495, "UC" 5644286,
               "A." 4008987, ".U" 3875698, "AC" 3812954, "CG" 3769734,
               "GU" 3742386, "G." 3532558, "C." 3352204, "GA" 340019,
               "CU" 315747, ".A" 288295, ".G" 243124, ".C" 142407,
               "GC" 72061, "U." 59768, "CA" 51845, "AU" 30890,
               "UG" 14227, "N." 353, "AN" 244, "CY" 178, "GN" 117,
               "CN" 99, ".S" 85, "NU" 75, "GS" 10, "SU" 9, "SA" 6,
               "CS" 4, ".N" 3, "NC" 2, "UN" 2, "NG" 1}]
       [fs (probs fs)]))

;;; (bg-freqs-probs 2 "/data2/Bio/RFAM"
;;;                 :fsps-fn cc-combins-freqs-probs
;;;                 :cols true :ftypes [".sto"] :sym? true :norm true :par 2)
(def +bg-freqs-probs-all-rfam-cols-sym+
     (let [fs {".." 165988848, "GG" 15598223, "AA" 14862080, "CC" 11859618,
               "UU" 8978659, "AG" 6959735, "UA" 6191385, "UC" 5960033,
               "A." 4297282, ".U" 3935466, "AC" 3864799, "CG" 3841795,
               "G." 3775682, "GU" 3756613, "C." 3494611, "N." 356,
               "AN" 244, "CY" 178, "GN" 118, "CN" 101, ".S" 85, "NU" 77,
               "GS" 10, "SU" 9, "SA" 6, "CS" 4}]
       [fs (probs fs)]))



(def
 ^{:doc
   "Memoized frequency and probability distributions over RFAM
    families.  The results returned are maps from an RFAM family
    identified by its string name (e.g., \"RF00050\") to pairs of the
    families frequency and probability distributions (as maps).  Keys
    are determined by n: 1 = mononucleotide, 2 = dinucleotide, etc.
    All nucleotide keys are strings.

    Ex:
    (sort (keys (bg-freqs-probs-rfam-by-family 1)))
    => ; A second or so
    (\"RF00050\" \"RF00057\" \"RF00059\" \"RF00114\" \"RF00140\" \"RF00162\"
     \"RF00167\" \"RF00168\" \"RF00174\" \"RF00234\" \"RF00379\" \"RF00380\"
     \"RF00442\" \"RF00504\" \"RF00506\" \"RF00514\" \"RF00515\" \"RF00521\"
     \"RF00522\" \"RF00555\" \"RF00556\" \"RF00557\" \"RF00558\" \"RF00559\"
     \"RF00634\" \"RF01051\" \"RF01054\" \"RF01055\" \"RF01057\" \"RF01070\"
     \"RF01385\" \"RF01402\" \"RF01482\" \"RF01510\" \"RF01689\" \"RF01692\"
     \"RF01693\" \"RF01694\" \"RF01704\" \"RF01725\" \"RF01727\" \"RF01739\"
     \"RF01745\" \"RF01767\" \"RF01769\" \"RF01793\" \"RF01826\" \"RF01831\")

     ((bg-freqs-probs-rfam-by-family 1) \"RF00050\")
     => ; Second time - cached and immediate
     [{\"C\" 4398, \"A\" 5034, \"G\" 6278, \"U\" 4212, \".\" 12344}
      {\".\" 0.3825698878075993, \"U\" 0.1305398871877518,
       \"G\" 0.1945701357466063, \"A\" 0.1560156201574413,
       \"C\" 0.1363044691006013}]

     (jensen-shannon
       (second ((bg-freqs-probs-rfam-by-family 2 :cols true) \"RF00050\"))
       (second +bg-freqs-probs-all-rfam-cols+))
     => 0.11288
   "

       :arglists '[n & {:keys [fsps-fn cols sym?]
                         :or {fsps-fn cc-freqs-probs cols false sym? false}}]}

     bg-freqs-probs-rfam-by-family

     (memoize
      (fn [n & {:keys [fsps-fn cols sym?]
                :or {fsps-fn cc-freqs-probs cols false sym? false}}]
        (if (and (= n 2)
                 (= fsps-fn cc-combins-freqs-probs)
                 cols)
          (io/with-in-reader "/data2/BioData/Archives/rfam-freqs-probs.clj"
            (read))
          (into
           {} (map (fn[file]
                     (let [nm (re-find #"^RF[0-9]+" (fs/basename file))
                           sqs (read-aln-seqs file :cols cols)
                           [fs ps] (take 2 (rest (seqs-freqs-probs
                                                  n sqs :fsps-fn fsps-fn
                                                  :nogaps false
                                                  :norm true :par 2)))]
                       [nm [fs ps]]))
                   (sort (fs/directory-files "/data2/Bio/RFAM" ".sto"))))))))


;;;(bg-freqs-probs-rfam-by-family 2)
;;;((bg-freqs-probs-rfam-by-family 2) "RF00050")
;;;
;;;
;;; This next one is quite expensive - several minutes and so we
;;; probably should secondarily cache this - that's where the DB
;;; should come in, but maybe just hack it to start with by writing it
;;; to a file...
;;;
;;;(bg-freqs-probs-rfam-by-family 2 :cols true :fsps-fn cc-combins-freqs-probs)

;;;(prn-dist
;;; (bg-freqs-probs-rfam-by-family 2 :cols true :fsps-fn cc-combins-freqs-probs)
;;; "/data2/BioData/Archives/rfam-freqs-probs.clj")

;;; (map #(do [% (jensen-shannon (second ((bg-freqs-probs-rfam-by-family
;;;                                        2 :cols true
;;;                                        :fsps-fn cc-combins-freqs-probs) %))
;;;                              (second +bg-freqs-probs-all-rfam-cols+))])
;;;      (sort (keys (bg-freqs-probs-rfam-by-family 1))))

;;; (def +ecr-rfam-js+
;;;      (for [sto (fs/directory-files "/data2/Bio/ECRibLeaders/MoStos" ".sto")
;;;            :let [nm (first (str/split #"\." (fs/basename sto)))
;;;                  stodist (second (bg-freqs-probs
;;;                                   2 [sto] :cols true
;;;                                   :fsps-fn cc-combins-freqs-probs
;;;                                   :par 2))]]
;;;        [nm (sort-by
;;;             second
;;;             (map #(do [% (jensen-shannon
;;;                           stodist
;;;                           (second ((bg-freqs-probs-rfam-by-family
;;;                                     2 :cols true
;;;                                     :fsps-fn cc-combins-freqs-probs) %)))])
;;;                  (sort (keys (bg-freqs-probs-rfam-by-family 1)))))]))

;;; (time
;;;  (bg-freqs-probs
;;;   2 (take 1 (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa"))
;;;   :par 4))


#_(time
   (io/with-out-writer "/data2/Bio/MetaG1/lHmax.clj"
     (doseq [mg-fa (take 20 (sort (fs/directory-files
                                   "/data2/Bio/MetaG1/FastaFiles" ".fa")))]
       (prn [(fs/basename mg-fa)
             (map #(let [[nm sq] %]
                     [nm (count sq) (-> sq count log4)
                      (for [l (range 1 20)]
                        [l (count (keys (reduce
                                         (fn[m [k v]]
                                           (if (> v 1)
                                             (assoc m k v) m))
                                         {} (freqn l sq))))])])
                  (read-seqs mg-fa :info :both))]))))


(comment

(def l-1 (probs 3 "CAAAAAGAAAUAGUAGGGCUUUUGAAAGUAACGGCCGCCCUGCAAAGAGCGGCCGUUCGCCUUACGUUUUUUA"))

(def l-2 (probs 2 "CAAAAAGAAAUAGUAGGGCUUUUGAAAGUAACGGCCGCCCUGCAAAGAGCGGCCGUUCGCCUUACGUUUUUUA"))

(for [x (keys l-1) a "AGUC"]
  (str x a))

(reduce (fn[m lmer]
          (let [l (dec (count lmer))
                x (subs lmer 1)
                y (subs lmer 0 l)
                z (subs lmer 1 l)]
            (if (and (l-1 x) (l-1 y) (l-2 z))
              (assoc m lmer (/ (* (l-1 x) (l-1 y)) (l-2 z)))
              m)))
        {} *1)

(jensen-shannon
 *1 (probs 4 "CAAAAAGAAAUAGUAGGGCUUUUGAAAGUAACGGCCGCCCUGCAAAGAGCGGCCGUUCGCCUUACGUUUUUUA"))

)


(defn expected-qdict [q-1 q-2 & {:keys [alpha] :or {alpha (alphabet :rna)}}]
  (reduce (fn[m lmer]
            (let [l (dec (count lmer))
                  x (subs lmer 1)
                  y (subs lmer 0 l)
                  z (subs lmer 1 l)]
              (if (and (q-1 x) (q-1 y) (q-2 z))
                (assoc m lmer (/ (* (q-1 x) (q-1 y)) (q-2 z)))
                m)))
          {} (for [x (keys q-1) a alpha] (str x a))))


(defn freq-xdict-dict
  [q sq]
  (let [ext-sq (str sq "X")
        residue-key (subs ext-sq (- (count sq) (- q 1)))
        q-xfreq-dict (freqn q ext-sq)
        q-freq-dict (dissoc q-xfreq-dict residue-key)]
    [q-xfreq-dict q-freq-dict]))

(defn q-1-dict
  ([q-xdict]
     (reduce (fn[m [k v]]
               (let [l (count k)
                     pre (subs k 0 (dec l))]
                 (assoc m pre (+ (get m pre 0) v))))
             {} q-xdict))
  ([q sq]
     (probs (dec q) sq)))

(defn q1-xdict-dict
  [q sq & {:keys [ffn] :or {ffn probs}}]
  (let [[q-xfreq-dict q-freq-dict] (freq-xdict-dict q sq)
        q-xpdf-dict (probs q-xfreq-dict)
        q-pdf-dict (probs q-freq-dict)]
    {:xfreq q-xfreq-dict :xpdf q-xpdf-dict
     :freq q-freq-dict :pdf q-pdf-dict}))



(defn reconstruct-dict
  [l sq]
  {:pre [(> l 2)]}
  (let [q (dec l)
        qmaps (q1-xdict-dict q sq)
        [q-xdict q-dict] (map qmaps [:xpdf :pdf])
        q-1dict (q-1-dict q-xdict)]
    (expected-qdict q-dict q-1dict)))


(defn max-qdict-entropy
  [q & {:keys [alpha] :or {alpha (alphabet :rna)}}]
  (let [base (count alpha)]
    (* q (log2 base))))

(defn informativity
  ([q sq]
     (- (max-qdict-entropy q) (entropy (probs q sq))))
  ([q-dict]
     (let [q (count (first (keys q-dict)))]
       (- (max-qdict-entropy q) (entropy q-dict)))))


(defn limit-entropy
  [q|q-dict sq|q-1dict &
   {:keys [alpha NA] :or {alpha (alphabet :rna) NA -1.0}}]
  {:pre [(or (and (integer? q|q-dict)
                  (or (string? sq|q-1dict) (coll? sq|q-1dict)))
             (and (map? q|q-dict) (map? sq|q-1dict)))]}

  (if (map? q|q-dict)
    (let [q-dict q|q-dict
          q-1dict sq|q-1dict]
      (/ (- (entropy q-dict) (entropy q-1dict))
         (log2 (count alpha))))

    (let [q q|q-dict
          sq sq|q-1dict
          lgcnt (log2 (count alpha))]
      (if (= q 1)
        (/ (entropy (probs 1 sq)) lgcnt)

        (let [qmaps (q1-xdict-dict q sq)
              [q-xdict q-dict] (map qmaps [:xpdf :pdf])
              q-1dict (q-1-dict q-xdict)]
          (if (< (count q-dict) (count q-1dict))
            (if (fn? NA) (NA q-dict q-1dict) NA)
            (/ (- (entropy q-dict) (entropy q-1dict)) lgcnt)))))))


(defn limit-informativity
  ([q sq]
     )
  ([q-dict]
     ))


(defn CREl [l sq & {:keys [limit] :or {limit 15}}]
  {:pre [(> l 2)]}
  (sum (fn[k]
         (catch-all (DX||Y
                     (probs k sq)
                     (reconstruct-dict k sq))))
       (range l (inc limit))))


(defn information-capacity
  [q sq & {:keys [cmpfn] :or {cmpfn jensen-shannon}}]
  (catch-all (cmpfn (probs q sq)
                    (reconstruct-dict q sq))))


(defn hybrid-dictionary
  [l sqs & {:keys [par] :or {par 1}}]
  {:pre [(or (string? sqs) (coll? sqs))]}
  (let [sqs (degap-seqs (if (coll? sqs) sqs (read-seqs sqs)))
        cnt (count sqs)
        dicts (pxmap #(probs l %) par sqs)
        Omega (apply set/union (map #(-> % keys set) dicts))
        hybrid (reducem (fn[k d] {k (get d k 0)})
                        #(merge-with + %1 %2)
                        Omega dicts)]
    (reduce (fn[m [k v]] (assoc m k (double (/ v cnt))))
            {} hybrid)))




#_(def *sto-hybrid-s4-0712*
     (let [l 7
           name-seqs (->> "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_0712.sto.S10regionprots.fna.cmsearch.csv" edu.bc.bio.gaisr.post-db-csv/get-entries
                          (keep (fn[[nm s e & _]]
                                  (when (fs/exists? (fs/join default-genome-fasta-dir (str nm ".fna")))
                                    (str nm "/" (if (> (Integer. s) (Integer. e)) (str e "-" s "/-1") (str s "-" e "/1"))))))
                          (pxmap gen-name-seq 10))
           names (map first name-seqs)
           sqs (map second name-seqs)
           sqs-hd (hybrid-dictionary l sqs :par 10)
           refsqs (->> (read-seqs "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_0712.sto") degap-seqs)
           sto-hd (hybrid-dictionary l refsqs :par 10)
           dicts (pxmap #(probs l %) 10 sqs)
           ;;sqsds (sort (pxmap #(DX||Y % sqs-hd) 10 dicts))
           stods (pxmap #(DX||Y % sto-hd) 10 dicts)
           nms-stods (sort-by second (map (fn[nm d] [nm d]) names stods))
           xs (range (count sqs))
           chart (incanter.charts/scatter-plot
                  xs (map second nms-stods)
                  :x-label "Sequence"
                  :y-label "RE to Hybrid"
                  :title (str "Sto to Hybrid, Resolution: " l))
           ;;chart (incanter.charts/add-points chart xs sqsds)
           ]
       (incanter.core/view chart)
       nms-stods))


#_(def *sto-hybrid-s4-0712-wctx*
     (let [l 7
           delta 5000
           names (->> "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_0712.sto.S10regionprots.fna.cmsearch.csv"
                      edu.bc.bio.gaisr.post-db-csv/get-entries
                      (keep (fn[[nm s e & _]]
                              (when (fs/exists? (fs/join default-genome-fasta-dir (str nm ".fna")))
                                (str nm "/" (if (> (Integer. s) (Integer. e)) (str e "-" s "/-1") (str s "-" e "/1")))))))
           name-seqs (->> names
                          (pxmap entry-parts 10)
                          (keep (fn[[nm [s e] sd]]
                                  (let [s (- s delta) s (if (<= s 0) 1 s) e (+ e delta)]
                                    (str nm "/" (if (> s e) (str e "-" s "/-1") (str s "-" e "/1"))))))
                          (pxmap gen-name-seq 10))
           sqs (map second name-seqs)
           sqs-hd (hybrid-dictionary l sqs :par 10)
           refsqs (->> (read-seqs "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_0712.sto" :info :name)
                       (pxmap entry-parts 10)
                       (pxmap (fn[[nm [s e] sd]]
                                (let [s (- s delta) s (if (<= s 0) 1 s) e (+ e delta)]
                                  (str/join "/" [nm (str s "-" e) sd]))) 10)
                       (pxmap gen-name-seq 10)
                       (map second))
           sto-hd (hybrid-dictionary l refsqs :par 10)
           dicts (pxmap #(probs l %) 10 sqs)
           ;;sqsds (sort (pxmap #(DX||Y % sqs-hd) 10 dicts))
           stods (pxmap #(DX||Y % sto-hd) 10 dicts)
           nms-stods (sort-by second (map (fn[nm d] [nm d]) names stods))
           xs (range (count sqs))
           chart (incanter.charts/scatter-plot
                  xs (map second nms-stods)
                  :x-label "Sequence"
                  :y-label "RE to Hybrid"
                  :title (str "Sto to Hybrid, Resolution: " l))
           ;;chart (incanter.charts/add-points chart xs sqsds)
           ]
       (incanter.core/view chart)
       nms-stods))

#_(let [candidates (->> (io/read-lines "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_072612a.ent")
                      (map entry-parts)
                      (map (fn[[nm [s e] sd]] [(str nm "/" s "-" e "/" sd) true]))
                      #_(take 15))
      s4-picked (into {} candidates)
      good (into {} (take 250 *sto-hybrid-s4-0712-wctx*))
      bad (into {} (drop 315 *sto-hybrid-s4-0712-wctx*))]
  [(count (reduce (fn[m [k v]] (if (get s4-picked k) (assoc m k v) m)) {} good)) (count s4-picked)])


#_(def sample-metag-cres
     (->> (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")
          first read-seqs
          (take 20)
          (map #(pxmap (fn[k] [k (CREl k %) (count %)]) 2 (range 3 16)))))

#_(doseq [cres (take 3 sample-metag-cres)]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)]
    (incanter.core/view
     (incanter.charts/scatter-plot
      ks cres
      :x-label "feature length"
      :y-label "Commulative Relative Entropy"
      :title (str "CRE(l, F/F^), len: " cnt
                  " Fmax: " (round (ln cnt)))))))

#_(doseq [cres foocres]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)]
    (incanter.core/view (incanter.charts/scatter-plot
                         ks cres
                         :x-label "feature length"
                         :y-label "Commulative Relative Entropy"
                         :title (str "CRE(l, F/F^), len: " cnt
                                     " Fmax: " (round (ln cnt)))))))

#_(def count-map
     (reduce (fn[m file]
               (reduce (fn[m sq]
                         (let [cnt (count sq)]
                           (assoc m cnt (inc (get m cnt 0)))))
                       m (read-seqs file)))
             {} (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")))


#_(def count-maps
       (map #(do [(fs/basename %)
                  (reduce (fn[m sq]
                            (let [cnt (count sq)]
                              (assoc m cnt (inc (get m cnt 0)))))
                          {} (read-seqs %))])
            (fs/directory-files "/data2/Bio/MetaG1/FastaFiles" ".fa")))

#_(sort-by
   first
   (map #(do [(first %) (->> % second (sort-by first >)
                             ((fn[sq]
                                [(take 5 sq) (drop (- (count sq) 5) sq)])))])
        count-maps))


#_(doseq [cres (map #(pxmap (fn[k] [k (CREl k %) (count %)]) 4 (range 3 16))
                  (->> "/home/jsa/Bio/STOfiles/S15stos/FreqDicts/S4_0712.sto.S10regionprots.fna.cmsearch.csv"
                       edu.bc.bio.gaisr.post-db-csv/get-entries
                       (keep (fn[[nm s e & _]]
                               (when (fs/exists? (fs/join default-genome-fasta-dir (str nm ".fna")))
                                 (str nm "/" (if (> (Integer. s) (Integer. e)) (str e "-" s "/-1") (str s "-" e "/1"))))))
                       gen-name-seq-pairs
                       (map second)
                       (take 7)))]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)
        Fmax (round (ln cnt))]
    (incanter.core/view
     (incanter.charts/scatter-plot
      ks cres
      :x-label "feature length"
      :y-label "Commulative Relative Entropy"
      :title (str "CRE(l, F/F^), len: " cnt
                  " Fmax: " Fmax)))))


#_(doseq [cres (map #(pxmap (fn[k] [k (CREl k %) (count %)]) 4 (range 3 16)) jsa)]
  (let [ks (map first cres)
        cnt (third (first cres))
        cres (map second cres)
        Fmax (round (ln cnt))]
    (incanter.core/save
     (incanter.charts/scatter-plot
      ks cres
      :x-label "feature length"
      :y-label "Commulative Relative Entropy"
      :title (str "CRE(l, F/F^), len: " cnt
                  " Fmax: " Fmax))
     (fs/join "/home/jsa/Bio/MetaGenome/Charts"
              (str "CREl-" cnt "-fmax-" Fmax ".png")))))
