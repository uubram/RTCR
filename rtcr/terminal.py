import time

def CUP(x,y):
    """Return terminal code to set cursor position to provided coordinates.
    
    :x: x-coordinate, 0 is left-most column
    :y: y-coordinate, 0 is top-most row
    """
    # ANSI has (1,1) as upper left corner
    return "\x1b[%s;%sH"%(y + 1, x + 1)

"""Return terminal code to set cursor to start of line."""
CUB = "\x1b[D"

def CLRSCR():
    return "\x1b[2J"

def EL(n = 0):
    return "\x1b[%sK"%n

def CHA(n):
    return "\x1b[%sG"%n

def CNL(n = 1):
    return "\x1b[%sE"%n

def CPL(n = 1):
    """Return terminal code to move cursor to beginning of line and up n
    lines.
    """
    return "\x1b[%sF"%n

def SU(n = 1):
    """Return terminal code to scroll up n lines."""
    return "\x1b[%sS"%n

class ProgressBar(object):
    def __init__(self, width = 40):
        self._width = width
        self.frac = 0.0
        self._start_time = time.time()

    def reset_start_time():
        self._start_time = time.time()

    def __str__(self):
        frac = self.frac
        assert 0 <= frac <= 1
        iwidth = self._width - 2

        if frac > 0:
            te = time.time() - self._start_time
            tr = (te / frac) * (1 - frac)
            trh = int(tr // 3600)
            trm = int((tr // 60) % 60)
            trs = int(tr % 60)
            s_remain = "(%sh:%sm:%ss remaining)"%(trh, trm, trs)
        else:
            s_remain = ""

        s_perc = "%.2f%%"%(frac * 100)

        if len(s_perc) + len(s_remain) > iwidth:
            if len(s_perc) > iwidth:
                s_inside = ""
            else:
                s_inside = s_perc
        else:
            s_inside = s_perc + s_remain

        # center the s_inside text
        n_blank = iwidth - len(s_inside)
        s_inside = (n_blank // 2)*" " + s_inside
        s_inside += " "*(iwidth - len(s_inside))

        n = int(frac * iwidth)
        s = "\x1b[30m\x1b[42m" + s_inside[:n] + "\x1b[0m" + s_inside[n:]
        return "[" + s + "]"
