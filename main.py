import schrodinger

func = schrodinger.WaveFunc(500, 0.002, 0.5, 10, -5, 25)
func.psi2_simulate_and_save(1, "psi2_k=0.5,V=0.gif", 0)
func.psi2_simulate_and_save(1, "psi2_k=0.5,V=gaus.gif", 1)

func.psi_simulate_and_save(1, "psi_k=0.5,V=0.gif", 0)
func.psi_simulate_and_save(1, "psi_k=0.5,V=gaus.gif", 1)

func = schrodinger.WaveFunc(500, 0.002, 0.75, 10, -5, 25)
func.psi2_simulate_and_save(1, "psi2_k=0.75,V=0.gif", 0)
func.psi2_simulate_and_save(1, "psi2_k=0.75,V=gaus.gif", 1)

func.psi_simulate_and_save(1, "psi_k=0.75,V=0.gif", 0)
func.psi_simulate_and_save(1, "psi_k=0.75,V=gaus.gif", 1)

func = schrodinger.WaveFunc(500, 0.002, 1, 10, -5, 25)
func.psi2_simulate_and_save(1, "psi2_k=1,V=0.gif", 0)
func.psi2_simulate_and_save(1, "psi2_k=1,V=gaus.gif", 1)

func.psi_simulate_and_save(1, "psi_k=1,V=0.gif", 0)
func.psi_simulate_and_save(1, "psi_k=1,V=gaus.gif", 1)

func = schrodinger.WaveFunc(500, 0.002, 1.5, 10, -5, 25)
func.psi2_simulate_and_save(1, "psi2_k=1.5,V=0.gif", 0)
func.psi2_simulate_and_save(1, "psi2_k=1.5,V=gaus.gif", 1)

func.psi_simulate_and_save(1, "psi_k=1,V=0.gif", 0)
func.psi_simulate_and_save(1, "psi_k=1,V=gaus.gif", 1)
