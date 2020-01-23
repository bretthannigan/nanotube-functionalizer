from Nanotube import Nanotube

n = Nanotube()
n.fromXYZ("<PATH_TO_CNT>.xyz")
n.functionalize("<PATH_TO_FUNC>.xyz", 0.1)
n.toXYZ("<PATH_TO_OUT>.xyz")