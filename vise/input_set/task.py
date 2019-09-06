from enum import unique, Enum


@unique
class Task(Enum):
    """ Supported tasks """
    single_point = "single_point"
    structure_opt = "structure_opt"
    phonon_disp = "phonon_disp"
    defect = "defect"
    band = "band"
    dos = "dos"
    dielectric = "dielectric"
    dielectric_function = "dielectric_function"
    gw_pre_calc1 = "gw_pre_calc1"
    gw_pre_calc2 = "gw_pre_calc2"

    def __str__(self):
        return self.name

    @classmethod
    def from_string(cls, s):
        for m in Task:
            if m.value == s:
                return m
        raise AttributeError("Task: " + str(s) + " is not proper.\n" +
                             "Supported Task:\n" + cls.name_list())

    @classmethod
    def name_list(cls):
        return ', '.join([e.value for e in cls])