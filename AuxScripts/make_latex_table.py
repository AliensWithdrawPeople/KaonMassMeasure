import numpy as np
import pandas as pd

energy = [501.0, 503.0, 505.0, 508.0, 508.5, 509.0, 509.5, 510.0, 510.5, 511.0, 511.5, 514.0, 517.0, 520.0, 525.0, 530.0]
systematics = [2.23, 0.86, 0.18, 0.02, 0.03, 0.68, 1.94, 0.75, 0.59, 0.02, 1.1, 0.02, 0.22, 1.06, 1.44, 0.71]
# End of setup.

energy = [f'${round(e, 1)}$' for e in energy]
systematics = [f'${round(e, 2)}$' for e in systematics]

df = pd.DataFrame({ r'$E_{beam}$, MeV': energy,
                    r'$\sigma_{syst}$, \%': systematics,
                    })  

with open('latex_EventCounting_systematics_table.txt', 'w+') as f:
    df.transpose().to_latex(f, 
                            index=True, 
                            header=False, 
                            column_format=f"{'|c' * (len(systematics) + 1)}|",
                            float_format="{:.1f}".format,
                            label='tab:EventCounting:Systematics', 
                            caption=r'Систематическая неопределённость подхода 2) процедуры подсчёта числа событий'
                            )
