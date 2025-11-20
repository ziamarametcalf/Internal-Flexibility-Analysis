python
from pymol import cmd

for i in range(1, 21):
    cmd.create(f"conf_{i}", "1anr", i, 1)
python end
