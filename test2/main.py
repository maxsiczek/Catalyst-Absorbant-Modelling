import CreateTrainingSet

c=CreateTrainingSet
c.sset_with_oxygen()
c.create_poscar_files('POTCAR_oxygen','structure_oxygen.pkl')
c.prepPoscar()

