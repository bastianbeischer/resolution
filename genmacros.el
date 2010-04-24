(defun copy-macro-files()
  (interactive)
  (cd "mac")
  (let ((angle 1.1))
    (while (< angle 3.0)
      (unless (= angle 1.0)
        (let ((new_angle (format "%.2f" angle)))
          (dolist (filename (directory-files "./" nil "perdaix_[0-9]\.[0-9]+_GeV_1\\.00_deg_msc_inhom\\.mac$" nil))
            (let ((new-file (replace-regexp-in-string "\\(.*\\)1\\.00\\(_deg_.*\\)" (format "\\1%s\\2" new_angle) filename)))
              (copy-file filename new-file t)
              (find-file new-file)
              (beginning-of-buffer)
              (if (re-search-forward "\\(moduleRot -\\)\\([0-9]+\\.[0-9]+\\)" nil t)
                  (replace-match (format "\\1%.2f" (/ angle 2.0))))
              (beginning-of-buffer)
              (if (re-search-forward "\\(moduleInternalRot \\)\\([0-9]+\\.[0-9]+\\)" nil t)
                  (replace-match (format "\\1%.2f" angle)))
              (beginning-of-buffer)
              (if (re-search-forward "\\(SetFileName .*\\)\\([0-9]+\\.[0-9]+\\)\\(_deg.*\\)" nil t)
                  (replace-match (format "\\1%.2f\\3" angle)))
              (save-buffer)
              (kill-buffer)))))
      (setq angle (+ angle 0.1)))))