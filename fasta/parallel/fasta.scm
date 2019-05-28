;; guile -e "(@ (fasta) fasta-cmd)" -s main.scm 25000000 > fasta.txt

(define-module (fasta)
  #:export (fasta-cmd fasta-main random-fasta)
  #:use-module (ice-9 iconv)
  #:use-module (ice-9 binary-ports)
  #:use-module (rnrs bytevectors)
  #:use-module (rnrs arithmetic fixnums)
  #:use-module (rnrs arithmetic flonums)
  ;; TODO parallel forms
  #:use-module (ice-9 threads))

(define (fasta-cmd args) (fasta-main (string->number (list-ref args 1))))

(define +alu+
  (string->bytevector
   (string-join
    '("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGG"
      "GAGGCCGAGGCGGGCGGATCACCTGAGGTCAGGAGTTCGAGA"
      "CCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAAT"
      "ACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCA"
      "GCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGG"
      "AGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCC"
      "AGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA")
    "")
   "iso-8859-1"))

(define (build-table t)
  (cons
   (u8-list->bytevector (map (compose char->integer car) t))
   (list->vector (map cdr t))))

(define IUB
  (build-table
   '([#\a . 0.27] [#\c . 0.12] [#\g . 0.12] [#\t . 0.27] [#\B . 0.02]
     [#\D . 0.02] [#\H . 0.02] [#\K . 0.02] [#\M . 0.02] [#\N . 0.02]
     [#\R . 0.02] [#\S . 0.02] [#\V . 0.02] [#\W . 0.02] [#\Y . 0.02])))

(define HOMOSAPIEN
  (build-table '([#\a . 0.3029549426680] [#\c . 0.1979883004921]
                 [#\g . 0.1975473066391] [#\t . 0.3015094502008])))

(define-values (IA IC IM LAST line-len) (values #f #f #f #f #f))

(define (init-fasta-vars)
  (set! IA 3877) (set! IC 29573) (set! IM 139968)
  (set! LAST 42) (set! line-len 60))

(define (random-next max)
  (set! LAST (modulo (+ IC (* LAST IA)) IM))
  (/ (* max (exact->inexact LAST)) (exact->inexact IM)))

(define (make-cumulative-table frequency-table)
  (define bs (car frequency-table))
  (define ps (cdr frequency-table))
  (define len (bytevector-length bs))
  (let loop ([i 0] [cum 0.0])
    (when (< i len)
      (let* ((this (vector-ref ps i))
             (new (+ this cum)))
        (vector-set! ps i new)
        (loop (1+ i) new)))))

(define (repeat-fasta s count port)
  (let* ((out port)
         (len (bytevector-length s))
         (len2 (+ len line-len))
         (s2 (make-bytevector len2)))
    (bytevector-copy! s 0 s2 0 len)
    (bytevector-copy! s 0 s2 len line-len)
    (let loop ([count count] [pos 0])
      (let* ((line (min line-len count))
             (count* (- count line)))
        (put-bytevector out s2 pos line)
        (newline out)
        (when (> count* 0)
          (let ((pos* (+ pos line)))
            (loop count* (if (>= pos* len) (- pos* len) pos*))))))))

(define (random-fasta genelist cnt port)
  (let ((out port)
        (ps (cdr genelist))
        (cs (car genelist)))
    (let loop ([count cnt])
      (let ((line (min line-len count))
            (buf (make-bytevector (1+ line-len))))
        (let inner ([pos 0])
          (let* ((r (random-next 1.0))
                 (i (let wh ([i 0]) (if (< (vector-ref ps i) r)
                                        (wh (1+ i))
                                        i))))
            (bytevector-u8-set! buf pos (bytevector-u8-ref cs i))
            (when (< (1+ pos) line)
              (inner (1+ pos)))))
        (bytevector-u8-set! buf line (char->integer #\newline))
        (put-bytevector out buf 0 (1+ line))
        (when (> (- count line) 0)
          (loop (- count line)))))))

(define (fasta-main n)
  (init-fasta-vars)
  (make-cumulative-table IUB)
  (make-cumulative-table HOMOSAPIEN)

  (let ((port (current-output-port)))
    (setvbuf port 'block (expt 2 24))
    (display ">ONE Homo sapiens alu\n" port)
    (repeat-fasta +alu+ (* n 2) port)
    (display ">TWO IUB ambiguity codes\n" port)
    (random-fasta IUB (* n 3) port)
    (display ">THREE Homo sapiens frequency\n" port)
    (random-fasta HOMOSAPIEN (* n 5) port)))
