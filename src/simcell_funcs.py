import os


def potcar(names,opsys):
    """
    names is a list of length m and each index is a string of the element name
    we will write to the POTCAR file in the order that is provided in the names
    list

    names = ['Type 1', 'Type 2',.., 'Type m']
    syntax of names should be captial letter followed by lowercase:

    He, Li, Au, etc
    """
    m = len(names)


    if opsys == "windows":
        beginning = 'type '
        file_name = ' potcar\{}_POTCAR'
    else:
        beginning = 'cat '
        file_name = ' potcar/{}_POTCAR'

    for i in range(m):
        ending = ' > POTCAR{}'.format(i+1)
        file_names = []
        for k in range(m):
            file_names.append(file_name.format(names[k]))
            middle = "".join(file_names)
            phrase = beginning+middle+ending

        os.system(phrase)


def formatter(alist):
    space1 = " "
    newline = '  '.join(str(num) for num in alist)+"\n"
    return(newline)
