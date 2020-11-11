import os
import subprocess


def removefile(what):
    try:
        os.remove(what)
    except OSError:
        pass


class OpsinConverter:

    def __init__(self, opsinjar, rm_tmp=True, wdir=os.path.abspath(os.getcwd()), silent=False):
        """
        a wrapper class for opsin

        :param opsinjar: absolute location for opsin.jar
        :param rm_tmp: if we should remove tmp files
        :param wdir: the working dir
        :param silent: if we should redirect stdout to null
        """
        self.opsinjar = opsinjar
        self.rm_tmp = rm_tmp
        self.wdir = wdir
        self.silent = silent

    def convert_name(self, name: str, output="smi"):
        tmpinfile = "opsin-" + str(hash(name)) + ".in"
        tmpoutfile = "opsin-" + str(hash(name)) + ".out"
        tmpinfile = "{}/{}".format(self.wdir, tmpinfile)
        tmpoutfile = "{}/{}".format(self.wdir, tmpoutfile)
        with open(tmpinfile, "w") as f:
            f.write(name)
        cmd = "java -jar {} -o{} {} {}".format(self.opsinjar, output, tmpinfile, tmpoutfile)
        if self.silent:
            proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        else:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        output = ""
        if os.path.isfile(tmpoutfile):
            with open(tmpoutfile, "r") as f:
                output = f.read().strip()
        if self.rm_tmp:
            removefile(tmpinfile)
            removefile(tmpoutfile)
        if len(output) == 0:
            return None
        return output

    def convert_names(self, names: [str], output="smi"):  # save i/o
        names = tuple(names)
        tmpinfile = "opsin-" + str(hash(names)) + ".in"
        tmpoutfile = "opsin-" + str(hash(names)) + ".out"
        tmpinfile = "{}/{}".format(self.wdir, tmpinfile)
        tmpoutfile = "{}/{}".format(self.wdir, tmpoutfile)
        with open(tmpinfile, "w") as f:
            for n in names:
                f.write(n + "\n")
        cmd = "java -jar {} -o{} {} {}".format(self.opsinjar, output, tmpinfile, tmpoutfile)
        if self.silent:
            proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        else:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        output = {}
        if os.path.isfile(tmpoutfile):
            with open(tmpoutfile, "r") as f:
                lines = f.readlines()
            for iline in range(len(lines)):
                line = lines[iline]
                if len(line.strip()) > 0:
                    output[names[iline]] = lines[iline].strip()
            # output = [l.strip() for l in lines if len(l.strip())>0]
        if self.rm_tmp:
            removefile(tmpinfile)
            removefile(tmpoutfile)
        return output
