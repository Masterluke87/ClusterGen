
def QMTemplates(prog = "ORCA", template="PBE-OPT",fp="./",coord="abc.xyz",charge=0,mult=1,procs=1):
    if prog=="ORCA":
        if template == "PBE-OPT":
            f = open(fp+"/pbe.inp","w")
            f.write(f"%pal nprocs {procs} end \n\n")
            f.write("!PBE def2-TZVP OPT UNO\n")
            f.write("!SLOWCONV\n\n")

            f.write(f"*xyzfile {charge} {mult}  {coord}")
            f.close()

            return True
        else:
            return False
    if prog=="GAUSSIAN":
        pass
    if prog=="PSI4":
        pass
    else:
        return False