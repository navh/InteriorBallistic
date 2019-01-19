from collections import namedtuple
from math import pi

BarrelResistance = namedtuple('BarrelResistance', 'br trav')
RecoilingPoint = namedtuple('RecoilingPoint', 'force time')
Propellant = namedtuple('Propellant', 'impetus flameTemp covolume mass density ratioOfSpecificHeats perforations lengthOfGrain diameterOfPerforation outsideDiameter burningRateList')
BurningRate = namedtuple('BurningRate', 'exponent coefficient pressure')

class main:
    def __init__(self,fn):
        self.fileName = fn

        self.readInAllData()

        self.f = open((self.fileName + '.out'), 'w+')
        self.printOutIdentifyingData()
        self.computeStoreConstantGroupings()
        self.run()

        self.f.close()

    def run(self):
        time = 0
        while(time <= float(self.time_to_stop)): #Make it end when the time ends
            if not ('Have individual propellants burned out?'):
                self.computeLinearBurningRates()
                self.computePropellantSurfaceAreas()
                self.computeMassFractionsBurnedByIntegration()
            if ('spacemeanpressure exceeded shot start pressure'):
                self.computeBreechPressure()
                self.computeBasePressure()
                self.interpolateForResistivePressure()
                self.computeProjectileAcceleration()
                self.computeProjectileVelocityByIntegration()
                self.computeProjectileDisplacementByIntegration()
            self.computeVolumeAvailableForPropellantGas()
            self.computeTemperatureOfPropellantGas()
            self.computeSpaceMeanPressure()
            self.checkForAndStoreMaxPressureAndAssociatedConditions()
            self.writeOutComputedResults()
            if not ('projectile has mooved'):
                if ('space mean pressure stopped increasing'):
                    break
            if ('projectile reached gun muzzle'):
                self.interpolateForConditionsAtMuzzle()
                self.writeConditionsAtMaximumPressureAndAtMuzzle()
                break
            time += float(self.time_step_sec)

        #TODO: page 02, pressure from the igniter Pa 185186.7

    def computeStoreConstantGroupings(self):
        self.calculateAreaOfTheBore()
        self.computeStapceMeanPressureAtTime0()

    def computeStapceMeanPressureAtTime0(self):
        # Pressure from the igniter, equation 42
        self.volume_of_unburnt_propellant = 0
        for propellant in self.propellants:
            self.volume_of_unburnt_propellant += (float(propellant.mass) / float(propellant.density))
        initial_volume = float(self.chamber_volume) - (float(self.covolume_of_igniter) * float(self.mass_of_igniter))
        self.igniter_pressure = float(self.impetus_of_igniter) * float(self.mass_of_igniter) / (initial_volume - float(self.volume_of_unburnt_propellant))
        self.f.write('pressure from the igniter Pa ' + str(self.igniter_pressure) + '\n')
        self.f.write('volume of unburnt propellant m^3 ' + str(self.volume_of_unburnt_propellant)+ '\n')
        self.f.write('initial chamber volume - covolume of ign m^3 ' + str(initial_volume) + '\n')

    def computeLinearBurningRates(self):
        # Using formula 32 (general case of formula 30)
        self.linear_burning_rates = []
        for propellant in self.propellants:
            linear_burning_rate = 0
            for br in propellant.burningRateList:
                pass
        pass

    def computePropellantSurfaceAreas(self):
        pass

    def computeMassFractionsBurnedByIntegration(self):
        pass

    def computeBreechPressure(self):
        pass

    def computeBasePressure(self):
        pass

    def interpolateForResistivePressure(self):
        pass

    def computeProjectileAcceleration(self):
        pass

    def computeProjectileVelocityByIntegration(self):
        pass

    def computeProjectileDisplacementByIntegration(self):
        pass

    def computeVolumeAvailableForPropellantGas(self):
        pass

    def computeTemperatureOfPropellantGas(self):
        pass

    def computeSpaceMeanPressure(self):
        pass

    def checkForAndStoreMaxPressureAndAssociatedConditions(self):
        pass

    def writeOutComputedResults(self):
        pass

    def interpolateForConditionsAtMuzzle(self):
        pass

    def writeConditionsAtMaximumPressureAndAtMuzzle(self):
        pass

    def readNextCaseOrStopProgram(self):
        pass

    def calculateAreaOfTheBore(self):
        grooveRadius = float(self.groove_diam) / 2
        landRadius = float(self.land_diameter) / 2
        grooveArea = pi * grooveRadius * grooveRadius
        landArea = pi * landRadius * landRadius
        sumOfRatio = (1 + float(self.groove_land_ratio))
        self.boreArea = ((grooveArea * float(self.groove_land_ratio)) + landArea) / sumOfRatio
        self.f.write('area of the bore m^2 ' + str(self.boreArea) + '\n')

    def readInAllData(self):
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

    def printOutIdentifyingData(self):
        self.f.write('the input file is ' + self.fileName + '.in' + '\n')
        self.f.write('the output file is ' + self.fileName + '.out' + '\n')
        self.f.write('using lagrange pressure gradient' + '\n')
        self.f.write(self.title)
        self.f.write('chamber volume in m^3 ' + self.chamber_volume + '\n')
        self.f.write('groove diam in m ' + self.groove_diam + '\n')
        self.f.write('land diameter in m ' + self.land_diameter + '\n')
        self.f.write('groove/land ratio ' + self.groove_land_ratio + '\n')
        self.f.write('twist in turns/caliber ' + self.twist_in_turns_caliber + '\n')
        self.f.write('travel in m ' + self.travel + '\n')
        self.f.write('gradient ' + self.gradient + '\n')
        self.f.write('' + '\n')
        self.f.write('projectile mass in kg ' + self.projectile_mass + '\n')
        self.f.write('switch to calculate if energy lost to air resistance ' + self.switch_to_calculate_energy_lost_to_air_resistance + '\n')
        self.f.write('fraction of work against bore used to heat tube ' + self.fraction_of_work_against_bore_to_heat_tube + '\n')
        self.f.write('gas pressure in front of projectile pa ' + self.gas_pressure_in_front_of_projectile + '\n')
        self.f.write('' + '\n')
        self.f.write('number of barrel resistance points (br,trav) ' + self.number_of_barrel_resistance_points + '\n')
        self.f.write('bore resistance Mpa     travel m' + '\n')
        for resistancePoint in self.barrel_resistance_points:
            self.f.write(' ' + resistancePoint.br + '\t\t\t\t\t\t\t' + resistancePoint.trav + '\n')
        self.f.write('' + '\n')
        self.f.write('mass of recoiling parts kg' + self.mass_of_recoiling_parts + '\n')
        self.f.write('number of recoil points (force,time) should be 2 ' + self.number_of_recoiling_parts + '\n')
        for part in self.recoiling_parts:
            self.f.write(' ' + part.force + '\t\t\t\t\t\t\t' + part.time + '\n')
        self.f.write('' + '\n')
        self.f.write('free convective heat transfer coefficient w/m^2-k ' + self.free_convective_heat_transfer_coefficient + '\n')
        self.f.write('chamber wall thickness m ' + self.chamber_wall_thickness + '\n')
        self.f.write('heat capacity of steel chamber wall j/kg-k ' + self.heat_capacity_of_steel_chamber_wall + '\n')
        self.f.write('initial temperature of chamber wall k ' + self.initial_temperature_of_chamber_wall + '\n')
        self.f.write('heat loss coefficient (should be 1) ' + self.heat_loss_coefficient + '\n')
        self.f.write('density of steel chamber wall kg/m^3 ' + self.density_of_steel_chamber_wall + '\n')
        self.f.write('' + '\n')
        self.f.write('impetus of igniter j/kg ' + self.impetus_of_igniter + '\n')
        self.f.write('covolume of igniter m^3/kg ' + self.covolume_of_igniter + '\n')
        self.f.write('adiabatic flame temperature k ' + self.adiabatic_flame_temperature + '\n')
        self.f.write('mass of igniter kg ' + self.mass_of_igniter + '\n')
        self.f.write('ratio of specific heats for igniter ' + self.ratio_of_specific_heats_for_igniter + '\n')
        self.f.write('' + '\n')
        self.f.write('number of propellants ' + self.number_of_propellants + '\n')
        self.f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            self.f.write('for propellant number ' + str(i) + '\n')
            self.f.write('impetus of the propellant j/kg ' + propellant.impetus + '\n')
            self.f.write('adiabatic flame temperature of propellant k ' + propellant.flameTemp + '\n')
            self.f.write('covolume of the propellant m^3/kg ' + propellant.covolume + '\n')
            self.f.write('mass of propellant kg ' + propellant.mass + '\n')
            self.f.write('density of propellant kg/m^3 ' + propellant.density + '\n')
            self.f.write('ratio of specific heats of propellant ' + propellant.ratioOfSpecificHeats + '\n')
            self.f.write('number of perforations of propellant grain ' + propellant.perforations + '\n')
            self.f.write('length of propellant grain m ' + propellant.lengthOfGrain + '\n')
            self.f.write('diameter of perforation of propellant grain m ' + propellant.diameterOfPerforation + '\n')
            self.f.write('outside diameter of propellant grain m ' + propellant.outsideDiameter + '\n')
        self.f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            self.f.write(' ' + str(len(propellant.burningRateList)) + ' burning rate points for propellant ' + str(i) + '\n')
            self.f.write('\n')
            self.f.write('exponent   coefficient  pressure' + '\n')
            self.f.write('  -         m/s-MPa^a      MPa' + '\n')
            for burnRate in propellant.burningRateList:
                self.f.write(' ' + burnRate.exponent + '\t\t'+ burnRate.coefficient + '\t\t' + burnRate.pressure + '\n')
        self.f.write('' + '\n')
        self.f.write('time step sec ' + self.time_step_sec + '\n')
        self.f.write('print step sec ' + self.print_step_sec + '\n')
        self.f.write('time to stop (if before projectile exit) sec ' + self.time_to_stop + '\n')

if __name__ == '__main__':
    main('19h')
