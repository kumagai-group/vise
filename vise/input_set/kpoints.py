# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp.sets import Kpoints


class ViseKpoints(Kpoints):

    def __str__(self):
        lines = [self.comment, str(self.num_kpts), self.style.name]
        style = self.style.name.lower()[0]
        if style == "l":
            lines.append(self.coord_type)
        for i in range(len(self.kpts)):
            # The following four lines are modified
            if len(self.kpts) == 1:
                lines.append(" ".join([str(x) for x in self.kpts[i]]))
            else:
                lines.append(" ".join([f"{x:20.17f}" for x in self.kpts[i]]))
            if style == "l":
                lines[-1] += " ! " + self.labels[i]
                if i % 2 == 1:
                    lines[-1] += "\n"
            elif self.num_kpts > 0:
                # The following four lines are modified
                if self.labels is not None:
                    lines[-1] += f"{self.kpts_weights[i]:4} {self.labels[i]}"
                else:
                    lines[-1] += f"{self.kpts_weights[i]:4}"

        # Print shifts for automatic kpoints types if not zero.
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            lines.append(" ".join([str(x) for x in self.kpts_shift]))
        return "\n".join(lines) + "\n"



