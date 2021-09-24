#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np

from pysimpp.fastpost import order_parameter # pylint: disable=no-name-in-module

class ClusterLog(object):
    def __init__(self, dirname, filename):
        self.dirname = dirname
        self.f = open(dirname + os.sep + filename, 'w')

    def get_file(self):
        return self.f

    def log(self, clusters):
        return

    def close(self):
        if not self.f is None:
            self.f.close()


class ClPropertiesLog(ClusterLog):
    def __init__(self, dirname, filename="properties.dat"):

        super(ClPropertiesLog, self).__init__(dirname, filename)

        header = "# " + " ".join(
            ("step", "nclusters", "nmolecules", "std", "sqrg", "std", "b",
             "std", "c", "std", "sqk", "std"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        # number of clusters
        nclusters = len(clusters)
        # number of molecules per cluster
        molecules = np.array([len(cl.molecules) for cl in clusters])
        # total number of molecules in the clusters
        nmolecules = np.sum(molecules)
        # radious of gyration per cluster
        sqrg = np.array([cl._sqrg for cl in clusters])
        # asphericity per cluster
        b = np.array([cl._b for cl in clusters])
        # acylindricity per cluster
        c = np.array([cl._c for cl in clusters])
        # shape anisotropy per clustesr
        sqk = np.array([cl._sqk for cl in clusters])

        self.f.write("%d %d %f %f %f %f %f %f %f %f %f %f\n" % # pylint: disable=bad-string-format-type
                     (step, len(clusters), molecules.mean(), molecules.std(),
                      sqrg.mean(), sqrg.std(), b.mean(), b.std(), c.mean(),
                      c.std(), sqk.mean(), sqk.std()))


class ClMolecularLog(ClusterLog):
    def __init__(self, dirname, filename="molecular.dat"):
        super(ClMolecularLog, self).__init__(dirname, filename)
        header = "# " + " ".join(("step", "nmolecules", "molsqee", "molsqrg",
                                  "molb", "molc", "molmsqk"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        # number of molecules per cluster
        molecules = np.array([len(cl.molecules) for cl in clusters])
        # total number of molecules in the clusters
        nmolecules = np.sum(molecules)
        # mean molecular square end-to-end of the molecules in the clusters
        msqee = np.array([cl._msqee * len(cl.molecules)
                          for cl in clusters]).sum() / nmolecules
        # mean molecular square radious of gyration of the molecules in the clusters
        msqrg = np.array([cl._msqrg * len(cl.molecules)
                          for cl in clusters]).sum() / nmolecules
        # mean molecular asphericity of the molecules in the clusters
        mb = np.array([cl._mb * len(cl.molecules)
                       for cl in clusters]).sum() / nmolecules
        # mean molecular acylindricity of the molecules in the clusters
        mc = np.array([cl._mc * len(cl.molecules)
                       for cl in clusters]).sum() / nmolecules
        # mean molecular shape anisotropy of the molecules in the clusters
        msqk = np.array([cl._msqk * len(cl.molecules)
                         for cl in clusters]).sum() / nmolecules

        self.f.write("%d %d %f %f %f %f %f\n" % # pylint: disable=bad-string-format-type
                     (step, nmolecules, msqee, msqrg, mb, mc, msqk))


class ClOrderLog(ClusterLog):
    def __init__(self, dirname, filename="order.dat"):
        super(ClOrderLog, self).__init__(dirname, filename)
        header = "# " + " ".join(
            ("step", "qloc", "std", "qlong", "std", "order[1]",
             "director[1][1]", "director[1][2]", "director[1][3]", "order[2]",
             "director[2][1]", "director[2][2]", "director[2][3]", "order[3]",
             "director[3][1]", "director[3][2]", "director[3][3]"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        # molecular local order paramters
        qloc = np.array([cl._qloc for cl in clusters])
        # order parameters
        qlong = np.array([cl._qlong for cl in clusters])
        # order parameter based on the primary cluster axis
        _v = np.array([cl._rgvec[0:3] for cl in clusters])
        _exclude = np.zeros(len(_v), dtype=np.bool)
        _q, _eigval, _eigvec, _ierr = order_parameter(_v, _exclude)

        self.f.write("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" % # pylint: disable=bad-string-format-type
                     (step, qloc.mean(), qloc.std(), qlong.mean(), qlong.std(),
                      _eigval[0], _eigvec[0], _eigvec[1], _eigvec[2],
                      _eigval[1], _eigvec[3], _eigvec[4], _eigvec[5],
                      _eigval[2], _eigvec[6], _eigvec[7], _eigvec[8]))


class ClDetailsLog(ClusterLog):
    def __init__(self, dirname, filename="details.dat"):
        super(ClDetailsLog, self).__init__(dirname, filename)
        header = "# " + " ".join(
            ("step", "id", "isinf", "nmols", "ldens", "ldens_std", "ldens_n", "ldens_L", "grval[1]",
             "grval[2]", "grval[3]", "bboxl[1]", "bbox[2]", "bbox[3]",
             "order[1]", "order[2]", "order[3]", "grvec[0][0]", "grvec[0][1]",
             "grvec[0][2]", "grvec[1][0]", "grvec[1][1]", "grvec[1][2]",
             "grvec[2][0]", "grvec[2][1]", "grvec[2][2]", "qvec[0][0]",
             "qvec[0][1]", "qvec[0][2]", "qvec[1][0]", "qvec[1][1]",
             "qvec[1][2]", "qvec[2][0]", "qvec[2][1]", "qvec[2][2]"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        for icl, cl in enumerate(clusters):
            nmols = len(cl.molecules)
            ldens = cl._ldensity
            bbox = cl._bbox
            grval = cl._rgval
            grvec = cl._rgvec
            qval = cl._qval
            qvec = cl._qvec
            inf = 1 if cl._infinit else 0
            self.f.write(
                "%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n"
                % (step, icl, inf, nmols, cl._ldensity[0], cl._ldensity[1], cl._ldensity[2], cl._ldensity[3],
                   grval[0], grval[1], grval[2], bbox[0], bbox[1], bbox[2],
                   qval[0], qval[1], qval[2], grvec[0], grvec[1], grvec[2],
                   grvec[3], grvec[4], grvec[5], grvec[6], grvec[7], grvec[8],
                   qvec[0], qvec[1], qvec[2], qvec[3], qvec[4], qvec[5],
                   qvec[6], qvec[7], qvec[8]))
