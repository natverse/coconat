# cosine_sim works

    Code
      cosine_sim(am)
    Output
                 1796817841 1797505019 1796818119 1827516355  818983130
      1796817841 1.00000000 0.09631851 0.00000000  0.0000000 0.03509041
      1797505019 0.09631851 1.00000000 0.00000000  0.0000000 0.09507157
      1796818119 0.00000000 0.00000000 1.00000000  0.0000000 0.07046832
      1827516355 0.00000000 0.00000000 0.00000000  1.0000000 0.20027384
      818983130  0.03509041 0.09507157 0.07046832  0.2002738 1.00000000

---

    Code
      pm <- prepare_cosine_matrix(cosine_sim(am, transpose = T))
