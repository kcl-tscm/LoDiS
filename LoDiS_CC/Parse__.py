    def parse_input(self):
        with open(self.filename, 'r') as f:
            data = f.readlines()
 
        self.n_atoms = int(data[0].split()[0])
        self.n_steps_total = int(len(data) / (self.n_atoms + 2))
        self.n_steps = self.n_steps_total // self.skip
         
        self.atom_list = [line.split()[0] for line in
                          data[2 : self.n_atoms + 2]]
 
        self.coordinates = sp.zeros((self.n_steps, self.n_atoms, 3))
        for step in range(self.n_steps):
            coords = sp.zeros((self.n_atoms, 3))
             
            i = step * self.skip * (self.n_atoms + 2)
            for j, line in enumerate(data[i + 2 : i + self.n_atoms + 2]):
                coords[j, :] = [float(value) for value in line.split()[1:]]
             
            self.coordinates[step] = coords