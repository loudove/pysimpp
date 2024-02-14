#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import abc

from pysimpp.fastpost import order_parameter  # pylint: disable=no-name-in-module


class ClusterLog(metaclass=abc.ABCMeta):
    ''' Prototype for a log keeper class. '''

    @abc.abstractmethod
    def __init__(self, dirname, filename):
        ''' Initialize a log file at the given direcory. '''
        self.dirname = dirname
        self.f = open(dirname + os.sep + filename, 'w')

    def get_file(self):
        ''' Return the file. '''
        return self.f

    @abc.abstractmethod
    def log(self, something):
        ''' Log something. '''
        pass

    def close(self):
        ''' Close the log file. '''
        if not self.f is None:
            self.f.close()


class ClPropertiesLog(ClusterLog):

    def __init__(self, dirname, filename="properties.dat"):
        ''' Initialize a log file for keeping clusters properties. '''

        super(ClPropertiesLog, self).__init__(dirname, filename)

        header = "# " + " ".join(
            ("step", "nclusters", "nagg", "std", "sqrg", "std", "b",
             "std", "c", "std", "sqk", "std"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        ''' Log cluster properties at the given step. '''

        if len(clusters) == 0:
            self.f.write("%d  0  %s" % (step, " 0.0"*10))
            return
        # number of clusters
        nclusters = len(clusters)
        # number of molecules per cluster
        molecules = np.array([len(cl.molecules) for cl in clusters])
        # total number of molecules in the clusters
        nmolecules = np.sum(molecules)
        # radious of gyration per cluster
        sqrg = np.array([cl.sqrg for cl in clusters])
        # asphericity per cluster
        b = np.array([cl.b for cl in clusters])
        # acylindricity per cluster
        c = np.array([cl.c for cl in clusters])
        # shape anisotropy per clustesr
        sqk = np.array([cl.sqk for cl in clusters])

        self.f.write("%d %d %f %f %f %f %f %f %f %f %f %f\n" %  # pylint: disable=bad-string-format-type
                     (step, len(clusters), molecules.mean(), molecules.std(),
                      sqrg.mean(), sqrg.std(), b.mean(), b.std(), c.mean(),
                      c.std(), sqk.mean(), sqk.std()))


class ClMolecularLog(ClusterLog):

    def __init__(self, dirname, moltypes, filename="molecular.dat"):
        ''' Initialize a log file for keeping cluster molecules properties. '''

        super(ClMolecularLog, self).__init__(dirname, filename)
        header = "# step nmolecules"
        for _mtype in moltypes:
            header += ' ' + \
                ' '.join(
                    [x+"_%s" % _mtype for x in ['n', 'sqee', 'sqrg', 'b', 'c', 'sqk']])
        self.f.write("%s\n" % header)

    def log(self, step, clusters, moltypes):
        ''' Log cluster molecules properties at the given step. '''

        # total number of molecules in the clusters
        nmolecules = np.sum(np.array([len(cl.molecules) for cl in clusters]))
        line = "%d %d " % (step, nmolecules)
        for _mtype in moltypes:
            # total number of molecules of type _mtype in the clusters
            nmol = np.sum(np.array([cl.nspecies[_mtype] for cl in clusters]))
            # mean molecular square end-to-end of the molecules in the clusters
            msqee = np.array([cl._msqee[_mtype] * cl.nspecies[_mtype]
                              for cl in clusters]).sum() / nmol if "_msqee" in vars(clusters[0]) else 0.0
            # mean molecular square radious of gyration of the molecules in the clusters
            msqrg = np.array([cl._msqrg[_mtype] * cl.nspecies[_mtype]
                              for cl in clusters]).sum() / nmol if "_msqrg" in vars(clusters[0]) else 0.0
            # mean molecular asphericity of the molecules in the clusters
            mb = np.array([cl._mb[_mtype] * cl.nspecies[_mtype]
                           for cl in clusters]).sum() / nmol if "_mb" in vars(clusters[0]) else 0.0
            # mean molecular acylindricity of the molecules in the clusters
            mc = np.array([cl._mc[_mtype] * cl.nspecies[_mtype]
                           for cl in clusters]).sum() / nmol if "_mc" in vars(clusters[0]) else 0.0
            # mean molecular shape anisotropy of the molecules in the clusters
            msqk = np.array([cl._msqk[_mtype] * cl.nspecies[_mtype]
                            for cl in clusters]).sum() / nmol if "_msqk" in vars(clusters[0]) else 0.0
            line += '%d %f %f %f %f %f ' % (
                nmol, msqee, msqrg, mb, mc, msqk)  # pylint: disable=bad-string-format-type
        self.f.write(line + '\n')


class ClOrderLog(ClusterLog):

    def __init__(self, dirname, filename="order.dat"):
        ''' Initialize a log file for keeping cluster order parameters. '''

        super(ClOrderLog, self).__init__(dirname, filename)
        header = "# " + " ".join(
            ("step", "qlocal", "std", "qlong", "std", "dirval[0]",
             "dirvec[0][0]", "dirvec[0][1]", "dirvec[0][2]", "dirval[1]",
             "dirvec[1][0]", "dirvec[1][1]", "dirvec[1][2]", "dirval[2]",
             "dirvec[2][0]", "dirvec[2][1]", "dirvec[2][2]"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        ''' Log cluster order parameters at the given step. '''

        if len(clusters) == 0:
            self.f.write("%d %s" % (step, " 0.0"*16))
            return
        # molecular local order paramters
        qlocal = np.array([cl.qlocal for cl in clusters])
        # order parameters
        qlong = np.array([cl.qlong for cl in clusters])
        # order parameter based on the primary cluster axis
        _v = np.array([cl.rgvec[0:3] for cl in clusters])
        _exclude = np.zeros(len(_v), dtype=np.bool_)
        _q, _eigval, _eigvec, _ierr = order_parameter(_v, _exclude)

        self.f.write("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" %  # pylint: disable=bad-string-format-type
                     (step, qlocal.mean(), qlocal.std(), qlong.mean(), qlong.std(),
                      _eigval[0], _eigvec[0], _eigvec[1], _eigvec[2],
                      _eigval[1], _eigvec[3], _eigvec[4], _eigvec[5],
                      _eigval[2], _eigvec[6], _eigvec[7], _eigvec[8]))


class ClDetailsLog(ClusterLog):

    def __init__(self, dirname, hasends, filename="details.dat"):
        ''' Initialize a log file for keeping system details. '''

        super(ClDetailsLog, self).__init__(dirname, filename)
        header = "# " + " ".join(
            ("step", "id", "isinf", "nagg", "ldens", "ldens_std", "ldens_n",
             "ldens_L", "bboxl[0]", "bbox[1]", "bbox[2]", "sval[0]", "sval[1]",
             "sval[2]", "svec[0][0]", "svec[0][1]", "svec[0][2]", "svec[1][0]",
             "svec[1][1]", "svec[1][2]", "svec[2][0]", "svec[2][1]",
             "svec[2][2]"))
        self.printq = hasends
        if hasends:
            header += " ".join(("qval[0]", "qval[1]", "qval[2]", "qvec[0][0]",
                                "qvec[0][1]", "qvec[0][2]", "qvec[1][0]", "qvec[1][1]",
                                "qvec[1][2]", "qvec[2][0]", "qvec[2][1]", "qvec[2][2]"))
        self.f.write("%s\n" % header)

    def log(self, step, clusters):
        ''' Log system details at the given step. '''

        for icl, cl in enumerate(clusters):
            nmols = len(cl.molecules)
            ldens = cl._ldensity
            bbox = cl.bbox
            grval = cl.rgval
            grvec = cl.rgvec
            inf = 1 if cl._infinit else 0
            line = "%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f" \
                % (step, icl, inf, nmols, cl._ldensity[0], cl._ldensity[1], cl._ldensity[2], cl._ldensity[3],
                   bbox[0], bbox[1], bbox[2],
                   grval[0], grval[1], grval[2], grvec[0], grvec[1], grvec[2],
                   grvec[3], grvec[4], grvec[5], grvec[6], grvec[7], grvec[8])

            if self.printq:
                qval = cl.qval
                qvec = cl.qvec
                line += " %f %f %f %f %f %f %f %f %f %f %f %f" \
                    % (qval[0], qval[1], qval[2], qvec[0], qvec[1], qvec[2],
                       qvec[3], qvec[4], qvec[5], qvec[6], qvec[7], qvec[8])
            self.f.write(line+"\n")
