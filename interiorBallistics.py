from collections import namedtuple
from math import pi

BarrelResistance = namedtuple('BarrelResistance', 'br trav')
RecoilingPoint = namedtuple('RecoilingPoint', 'force time')
Propellant = namedtuple('Propellant', 'impetus flameTemp covolume mass density ratioOfSpecificHeats perforations lengthOfGrain diameterOfPerforation outsideDiameter burningRateList')
BurningRate = namedtuple('BurningRate', 'exponent coefficient pressure')

class main:
    def __init__(self,fn):
        self.fileName = fn

        self.parseFile()
        self.calculateAreaOfTheBore()
        self.writeParseInfo()
        #TODO: page 02, pressure from the igniter Pa 185186.7

    def calculateAreaOfTheBore(self):
        grooveRadius = float(self.groove_diam) / 2
        landRadius = float(self.land_diameter) / 2
        grooveArea = pi * grooveRadius * grooveRadius
        landArea = pi * landRadius * landRadius
        sumOfRatio = (1 + float(self.groove_land_ratio))
        self.boreArea = ((grooveArea * float(self.groove_land_ratio)) + landArea) / sumOfRatio

    def parseFile(self):
        f = open((self.fileName + '.in'), "r")
        self.title = f.readline().replace('"', '')

        line = f.readline().split()
        self.chamber_volume = line[0]
        self.groove_diam = line[1]
        self.land_diameter = line[2]
        self.groove_land_ratio = line[3]
        self.twist_in_turns_caliber = line[4]
        self.travel = line[5]
        self.gradient = line[6]

        line = f.readline().split()
        self.projectile_mass = line[0]
        self.switch_to_calculate_energy_lost_to_air_resistance = line[1]
        self.fraction_of_work_against_bore_to_heat_tube = line[2]
        self.gas_pressure_in_front_of_projectile = line[3]

        line = f.readline().split()
        self.number_of_barrel_resistance_points = line[0]
        self.barrel_resistance_points = []
        for _ in range(int(self.number_of_barrel_resistance_points)):
            line = f.readline().split()
            self.barrel_resistance_points.append(BarrelResistance(line[0],line[1]))

        line = f.readline().split()
        self.mass_of_recoiling_parts = line[0]
        self.number_of_recoiling_parts = line[1]

        self.recoiling_parts = []
        for _ in range(int(self.number_of_recoiling_parts)):
            line = f.readline().split()
            self.recoiling_parts.append(RecoilingPoint(line[0],line[1]))

        line = f.readline().split()
        self.free_convective_heat_transfer_coefficient = line[0]
        self.chamber_wall_thickness = line[1]
        self.heat_capacity_of_steel_chamber_wall = line[2]
        self.initial_temperature_of_chamber_wall = line[3]
        self.heat_loss_coefficient = line[4]
        self.density_of_steel_chamber_wall = line[5]

        line = f.readline().split()
        self.impetus_of_igniter = line[0]
        self.covolume_of_igniter = line[1]
        self.adiabatic_flame_temperature = line[2]
        self.mass_of_igniter = line[3]
        self.ratio_of_specific_heats_for_igniter = line[4]

        line = f.readline().split()
        self.number_of_propellants = line[0]

        self.propellants = []
        for _ in range(int(self.number_of_propellants)):
            line = f.readline().split()
            self.propellants.append(Propellant(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],[]))

        for propellant in self.propellants:
            numberOfBurningRates = int(f.readline())
            for _ in range(numberOfBurningRates):
                line = f.readline().split()
                propellant.burningRateList.append(BurningRate(line[0],line[1],line[2]))

        line = f.readline().split()
        self.time_step_sec = line[0]
        self.print_step_sec = line[1]
        self.time_to_stop = line[2]

        f.close()

    def writeParseInfo(self):
        f = open((self.fileName + '.out'), 'w+')
        f.write('the input file is ' + self.fileName + '.in' + '\n')
        f.write('the output file is ' + self.fileName + '.out' + '\n')
        f.write('using lagrange pressure gradient' + '\n')
        f.write(self.title)
        f.write('chamber volume in m^3 ' + self.chamber_volume + '\n')
        f.write('groove diam in m ' + self.groove_diam + '\n')
        f.write('land diameter in m ' + self.land_diameter + '\n')
        f.write('groove/land ratio ' + self.groove_land_ratio + '\n')
        f.write('twist in turns/caliber ' + self.twist_in_turns_caliber + '\n')
        f.write('travel in m ' + self.travel + '\n')
        f.write('gradient ' + self.gradient + '\n')
        f.write('' + '\n')
        f.write('projectile mass in kg ' + self.projectile_mass + '\n')
        f.write('switch to calculate if energy lost to air resistance ' + self.switch_to_calculate_energy_lost_to_air_resistance + '\n')
        f.write('fraction of work against bore used to heat tube ' + self.fraction_of_work_against_bore_to_heat_tube + '\n')
        f.write('gas pressure in front of projectile pa ' + self.gas_pressure_in_front_of_projectile + '\n')
        f.write('' + '\n')
        f.write('number of barrel resistance points (br,trav) ' + self.number_of_barrel_resistance_points + '\n')
        f.write('bore resistance Mpa     travel m' + '\n')
        for resistancePoint in self.barrel_resistance_points:
            f.write(' ' + resistancePoint.br + '\t\t\t\t\t\t\t' + resistancePoint.trav + '\n')
        f.write('' + '\n')
        f.write('mass of recoiling parts kg' + self.mass_of_recoiling_parts + '\n')
        f.write('number of recoil points (force,time) should be 2 ' + self.number_of_recoiling_parts + '\n')
        for part in self.recoiling_parts:
            f.write(' ' + part.force + '\t\t\t\t\t\t\t' + part.time + '\n')
        f.write('' + '\n')
        f.write('free convective heat transfer coefficient w/m^2-k ' + self.free_convective_heat_transfer_coefficient + '\n')
        f.write('chamber wall thickness m ' + self.chamber_wall_thickness + '\n')
        f.write('heat capacity of steel chamber wall j/kg-k ' + self.heat_capacity_of_steel_chamber_wall + '\n')
        f.write('initial temperature of chamber wall k ' + self.initial_temperature_of_chamber_wall + '\n')
        f.write('heat loss coefficient (should be 1) ' + self.heat_loss_coefficient + '\n')
        f.write('density of steel chamber wall kg/m^3 ' + self.density_of_steel_chamber_wall + '\n')
        f.write('' + '\n')
        f.write('impetus of igniter j/kg ' + self.impetus_of_igniter + '\n')
        f.write('covolume of igniter m^3/kg ' + self.covolume_of_igniter + '\n')
        f.write('adiabatic flame temperature k ' + self.adiabatic_flame_temperature + '\n')
        f.write('mass of igniter kg ' + self.mass_of_igniter + '\n')
        f.write('ratio of specific heats for igniter ' + self.ratio_of_specific_heats_for_igniter + '\n')
        f.write('' + '\n')
        f.write('number of propellants ' + self.number_of_propellants + '\n')
        f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            f.write('for propellant number ' + str(i) + '\n')
            f.write('impetus of the propellant j/kg ' + propellant.impetus + '\n')
            f.write('adiabatic flame temperature of propellant k ' + propellant.flameTemp + '\n')
            f.write('covolume of the propellant m^3/kg ' + propellant.covolume + '\n')
            f.write('mass of propellant kg ' + propellant.mass + '\n')
            f.write('density of propellant kg/m^3 ' + propellant.density + '\n')
            f.write('ratio of specific heats of propellant ' + propellant.ratioOfSpecificHeats + '\n')
            f.write('number of perforations of propellant grain ' + propellant.perforations + '\n')
            f.write('length of propellant grain m ' + propellant.lengthOfGrain + '\n')
            f.write('diameter of perforation of propellant grain m ' + propellant.diameterOfPerforation + '\n')
            f.write('outside diameter of propellant grain m ' + propellant.outsideDiameter + '\n')
        f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            f.write(' ' + str(len(propellant.burningRateList)) + ' burning rate points for propellant ' + str(i) + '\n')
            f.write('\n')
            f.write('exponent   coefficient  pressure' + '\n')
            f.write('  -         m/s-MPa^a      MPa' + '\n')
            for burnRate in propellant.burningRateList:
                f.write(' ' + burnRate.exponent + '\t\t'+ burnRate.coefficient + '\t\t' + burnRate.pressure + '\n')
        f.write('' + '\n')
        f.write('time step sec ' + self.time_step_sec + '\n')
        f.write('print step sec ' + self.print_step_sec + '\n')
        f.write('time to stop (if before projectile exit) sec ' + self.time_to_stop + '\n')
        f.write('area of the bore m^2 ' + str(self.boreArea) + '\n')


        f.close()

if __name__ == '__main__':
    main('19h')
