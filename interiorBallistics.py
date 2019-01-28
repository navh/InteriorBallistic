from collections import namedtuple
from math import pi
from scipy.constants import g as gravity

PIDDUCK_KENT_CONSTANT = 1

BarrelResistance = namedtuple('BarrelResistance', 'br trav')
RecoilingPoint = namedtuple('RecoilingPoint', 'force time')
BurningRate = namedtuple('BurningRate', 'exponent coefficient pressure')


class Propellant:
    def __init__(self, imp, fla, cov,mass,dens,ratSpecHeat,perfs,lenGra,diaPerf,outDia):
        self.impetus = float(imp)
        self.flameTemp  = float(fla)
        self.covolume = float(cov)
        self.mass = float(mass)
        self.density  = float(dens)
        self.ratioOfSpecificHeats = float(ratSpecHeat)
        self.perforations = float(perfs)
        self.lengthOfGrain  = float(lenGra)
        self.diameterOfPerforation = float(diaPerf)
        self.outsideDiameter = float(outDia)
        self.burningRateList = []
        self.surfaceArea = None
        self.linearBurnRate = None
        self.massFractionBurningRate = None
        self.massFraction = None

    def initializePropellant(self):
        self.computeBurningRate()
        self.computeInitialSurfaceArea()
        self.computeMassFraction(0)


    def computeInitialSurfaceArea(self):
        self.surfaceArea = pi * (((self.outsideDiameter) + (self.perforations) * (self.diameterOfPerforation)) * (
            (self.lengthOfGrain)) + (((self.outsideDiameter) ** 2 - (
                (self.perforations) * ((self.diameterOfPerforation))) ** 2) / 2))

    def computeSurfaceArea(self,time):
        ui = 2 * (self.linearBurnRate) * time # Big guess for ui
        self.surfaceArea = pi * (((self.outsideDiameter) - ui + (self.perforations) * (
                    (self.diameterOfPerforation) + ui)) * (((self.lengthOfGrain) - ui)) + ((((
            self.outsideDiameter) - ui) ** 2 - ((self.perforations) * (
                    ((self.diameterOfPerforation) + ui) ** 2))) / 2))

    def computeMassFraction(self,time):
        ui = 2 * (self.linearBurnRate) * time  # Big guess for ui
        vgi = (pi / 4) * ((self.outsideDiameter) ** 2 - ((self.perforations) * ((self.diameterOfPerforation) ** 2))) * self.lengthOfGrain #Initial volume using geometry instead of density and mass
        volumeOfPartiallyBurntPropellant = (pi / 4) * (((self.outsideDiameter - ui)**2 - (self.perforations * (self.diameterOfPerforation + ui)**2)) * (self.lengthOfGrain - ui))
        self.massFraction = (1 - (volumeOfPartiallyBurntPropellant / vgi))
        self.massFractionBurningRate = 1 / (vgi * self.surfaceArea * self.linearBurnRate)

    def computeBurningRate(self):
        self.linearBurnRate = 0
        for br in self.burningRateList:
            self.linearBurnRate += (((br.coefficient) * (br.pressure)) ** (br.exponent))  # TODO: Check to see if summing these makes sense, we did to see if our example with only 1 burn rate / propellant would work

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
        while(self.t <= (self.time_to_stop)):
            if not ('Have individual propellants burned out?'):
                self.computeLinearBurningRates()
                self.computePropellantSurfaceAreas()
                self.computeMassFractionsBurnedByIntegration()
            if (self.space_mean_pressure > self.retardingPressureByMeter(0)):
                self.computeBasePressure()
                self.computeBreechPressure()
                self.interpolateForResistivePressure()
                self.computeProjectileAcceleration()
                self.computeProjectileVelocityByIntegration()
                self.computeProjectileDisplacementByIntegration()
            self.computeVolumeAvailableForPropellantGas()
            self.computeTemperatureOfPropellantGas()
            self.computeSpaceMeanPressure()
            self.checkForAndStoreMaxPressureAndAssociatedConditions()
            self.writeOutComputedResults()
            if not (self.displacement_of_projectile > 0.0):
                if ('space mean pressure stopped increasing'):
                    break
            if (self.displacement_of_projectile > self.travel):
                self.interpolateForConditionsAtMuzzle()
                self.writeConditionsAtMaximumPressureAndAtMuzzle()
                break
            self.t += self.time_step_sec

    def computeStoreConstantGroupings(self):
        self.t = self.time_step_sec
        self.calculateAreaOfTheBore()
        self.computeStapceMeanPressureAtTime0()
        for p in self.propellants:
            p.initializePropellant()
        self.velocity_of_projectile = 0.0
        self.acceleration_of_projectile = 0.0
        self.displacement_of_projectile = 0.0
        self.f.write('  time           acc            vel            dis            mpress         pbase          pbrch\n')
        self.f.write('   s              m/s^2          m/s            m              Pa             Pa             Pa \n')


    def initialPropellantSurfaceArea(self):
        for p in self.propellants:
            p.computeSurfaceArea(0) #TODO: Make sure that this is the correct way to do this

    def computeStapceMeanPressureAtTime0(self):
        # Pressure from the igniter, equation 42
        self.volume_of_unburnt_propellant = 0
        for propellant in self.propellants:
            self.volume_of_unburnt_propellant += ((propellant.mass) / (propellant.density))
        initial_volume = (self.chamber_volume) - ((self.covolume_of_igniter) * (self.mass_of_igniter))
        self.igniter_pressure = (self.impetus_of_igniter) * (self.mass_of_igniter) / (initial_volume - (self.volume_of_unburnt_propellant))
        self.space_mean_pressure = self.igniter_pressure
        self.f.write('pressure from the igniter Pa ' + str(self.igniter_pressure) + '\n')
        self.f.write('volume of unburnt propellant m^3 ' + str(self.volume_of_unburnt_propellant)+ '\n')
        self.f.write('initial chamber volume - covolume of ign m^3 ' + str(initial_volume) + '\n')

    def computeLinearBurningRates(self):
        # Using formula 32 (general case of formula 30)
        for propellant in self.propellants:
            propellant.computeBurningRate()

    def computePropellantSurfaceAreas(self):
        for p in self.propellants:
            p.computeSurfaceArea(self.t)

    def computeMassFractionsBurnedByIntegration(self):
        for p in self.propellants:
            p.computeMassFraction(self.t)

    def computeBasePressure(self):
        sumOfMasses = 0
        for p in self.propellants:
            sumOfMasses = p.mass

        self.base_pressure = self.space_mean_pressure / (
                    1 + (sumOfMasses / (self.projectile_mass * PIDDUCK_KENT_CONSTANT)))

    def computeBreechPressure(self):
        """
        Equation 28, it got ugly again so we did the A B C trick again.

        :return:
        """

        AA = 0
        AB = 0
        CA = 0
        for p in self.propellants:
            AA += (p.mass * p.ratioOfSpecificHeats)
            AB += (p.mass)
            CA += (p.mass / self.projectile_mass)

        A = AA / AB

        B = 1 / (A - 1)

        C = ((2 * B + 3) / PIDDUCK_KENT_CONSTANT) + ((2 * (B + 1)) / (CA))

        self.breech_pressure = self.base_pressure / ((1 - (1/C)) ** (-B - 1))

    def interpolateForResistivePressure(self):
        self.resistive_pressure = self.retardingPressureByMeter(self.displacement_of_projectile)

    def computeProjectileAcceleration(self):
        self.acceleration_of_projectile = (self.boreArea * gravity * (self.base_pressure - self.gas_pressure_in_front_of_projectile - self.resistive_pressure)) / self.projectile_mass #F = MA BAYBEEE

    def computeProjectileVelocityByIntegration(self):
        self.velocity_of_projectile += (self.acceleration_of_projectile * self.t - self.acceleration_of_projectile * (self.t - self.time_step_sec))

    def computeProjectileDisplacementByIntegration(self):
        self.displacement_of_projectile += (self.velocity_of_projectile * self.t - self.velocity_of_projectile * (self.t - self.time_step_sec))

    def computeVolumeAvailableForPropellantGas(self):
        volume_occupied_by_unburned_solid_propellant = 0
        volume_occupied_by_gas_molecules = 0
        for p in self.propellants:
            volume_occupied_by_unburned_solid_propellant += (p.mass / p.density) * (1 - p.massFraction)
            volume_occupied_by_gas_molecules += (p.mass * p.massFraction * p.covolume)
        self.volume_available_for_propellant_gas = self.chamber_volume + self.boreArea * self.displacement_of_projectile - volume_occupied_by_gas_molecules - volume_occupied_by_unburned_solid_propellant

    def computeTemperatureOfPropellantGas(self):
        """
        From Equation (19)

        Our documentation provides a very long formula. For readability it has been broken down into 7 sub-formulae A-G

        :return:
        """
        A = 0
        for p in self.propellants:
            A += (p.impetus * p.mass * p.massFraction) / (p.ratioOfSpecificHeats - 1)

        B = (self.impetus_of_igniter * self.mass_of_igniter)/ (self.ratio_of_specific_heats_for_igniter - 1)

        CSum = 0
        for p in self.propellants:
            CSum = p.mass
        C = (self.velocity_of_projectile ** 2 / (gravity * 2)) * (self.projectile_mass + CSum / PIDDUCK_KENT_CONSTANT)

        D = - self.boreArea * self.retardingPressureByMeter(self.displacement_of_projectile) * self.displacement_of_projectile #TODO: We used a random integral here and are just guessing that there is no need to evaluate it

        E = 0 #TODO: Heat is usually insignificant, the formula for its calculation is complex (formula 17), it may be added later

        F = 0
        for p in self.propellants:
            F += (p.impetus * p.mass * p.massFraction) / ((p.ratioOfSpecificHeats - 1) * p.flameTemp)

        G = (self.impetus_of_igniter * self.mass_of_igniter) / ((self.ratio_of_specific_heats_for_igniter - 1) * self.adiabatic_flame_temperature)

        self.temperature_of_propellant_gas = (A + B - C - D - E) / (F + G)

    def computeSpaceMeanPressure(self):
        """
        Equation 26
        :return:
        """

        A = 0
        for p in self.propellants:
            A += (p.impetus * p.mass * p.massFraction) / (p.flameTemp)

        B = self.impetus_of_igniter * self.mass_of_igniter / self.adiabatic_flame_temperature

        self.space_mean_pressure = (self.temperature_of_propellant_gas / self.volume_available_for_propellant_gas) * (A + B)

    def checkForAndStoreMaxPressureAndAssociatedConditions(self):
        pass #TODO: Instructions unclear, pass intentionally left in place

    def writeOutComputedResults(self):
        self.f.write("%08.8E %08.8E %08.8E %08.8E %08.8E %08.8E %08.8E \n" % (self.t, self.acceleration_of_projectile, self.velocity_of_projectile, self.displacement_of_projectile, self.space_mean_pressure,self.space_mean_pressure,self.space_mean_pressure))


    def interpolateForConditionsAtMuzzle(self):
        pass

    def writeConditionsAtMaximumPressureAndAtMuzzle(self):
        pass

    def readNextCaseOrStopProgram(self):
        pass

    def retardingPressureByMeter(self,travelInMeters):
        pressure = None
        distance = -1
        for barrelResistancePoint in self.barrel_resistance_points:
            if (barrelResistancePoint.trav) <= travelInMeters:
                if (barrelResistancePoint.trav) > (distance):
                    distance = (barrelResistancePoint.trav)
                    pressure = (barrelResistancePoint.br)
        return pressure


    def calculateAreaOfTheBore(self):
        grooveRadius = (self.groove_diam) / 2
        landRadius = (self.land_diameter) / 2
        grooveArea = pi * grooveRadius * grooveRadius
        landArea = pi * landRadius * landRadius
        sumOfRatio = (1 + (self.groove_land_ratio))
        self.boreArea = ((grooveArea * (self.groove_land_ratio)) + landArea) / sumOfRatio
        self.f.write('area of the bore m^2 ' + str(self.boreArea) + '\n')

    def readInAllData(self):
        f = open((self.fileName + '.in'), "r")
        self.title = f.readline().replace('"', '')

        line = f.readline().split()
        self.chamber_volume = float(line[0])
        self.groove_diam = float(line[1])
        self.land_diameter = float(line[2])
        self.groove_land_ratio = float(line[3])
        self.twist_in_turns_caliber = float(line[4])
        self.travel = float(line[5])
        self.gradient = float(line[6])

        line = f.readline().split()
        self.projectile_mass = float(line[0])
        self.switch_to_calculate_energy_lost_to_air_resistance = float(line[1])
        self.fraction_of_work_against_bore_to_heat_tube = float(line[2])
        self.gas_pressure_in_front_of_projectile = float(line[3])

        line = f.readline().split()
        self.number_of_barrel_resistance_points = float(line[0])
        self.barrel_resistance_points = []
        for _ in range(int(self.number_of_barrel_resistance_points)):
            line = f.readline().split()
            self.barrel_resistance_points.append(BarrelResistance(float(line[0]),float(line[1])))

        line = f.readline().split()
        self.mass_of_recoiling_parts = float(line[0])
        self.number_of_recoiling_parts = float(line[1])

        self.recoiling_parts = []
        for _ in range(int(self.number_of_recoiling_parts)):
            line = f.readline().split()
            self.recoiling_parts.append(RecoilingPoint(float(line[0]),float(line[1])))

        line = f.readline().split()
        self.free_convective_heat_transfer_coefficient = float(line[0])
        self.chamber_wall_thickness = float(line[1])
        self.heat_capacity_of_steel_chamber_wall = float(line[2])
        self.initial_temperature_of_chamber_wall = float(line[3])
        self.heat_loss_coefficient = float(line[4])
        self.density_of_steel_chamber_wall = float(line[5])

        line = f.readline().split()
        self.impetus_of_igniter = float(line[0])
        self.covolume_of_igniter = float(line[1])
        self.adiabatic_flame_temperature = float(line[2])
        self.mass_of_igniter = float(line[3])
        self.ratio_of_specific_heats_for_igniter = float(line[4])

        line = f.readline().split()
        self.number_of_propellants = float(line[0])

        self.propellants = []
        for _ in range(int(self.number_of_propellants)):
            line = f.readline().split()
            self.propellants.append(Propellant(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]))

        for propellant in self.propellants:
            numberOfBurningRates = int(f.readline())
            for _ in range(numberOfBurningRates):
                line = f.readline().split()
                propellant.burningRateList.append(BurningRate(float(line[0]),float(line[1]),float(line[2])))

        line = f.readline().split()
        self.time_step_sec = float(line[0])
        self.print_step_sec = float(line[1])
        self.time_to_stop = float(line[2])

        f.close()

    def printOutIdentifyingData(self):
        self.f.write('the input file is ' + str(self.fileName) + '.in' + '\n')
        self.f.write('the output file is ' + str(self.fileName) + '.out' + '\n')
        self.f.write('using lagrange pressure gradient' + '\n')
        self.f.write(self.title)
        self.f.write('chamber volume in m^3 ' + str(self.chamber_volume) + '\n')
        self.f.write('groove diam in m ' + str(self.groove_diam) + '\n')
        self.f.write('land diameter in m ' + str(self.land_diameter) + '\n')
        self.f.write('groove/land ratio ' + str(self.groove_land_ratio) + '\n')
        self.f.write('twist in turns/caliber ' + str(self.twist_in_turns_caliber) + '\n')
        self.f.write('travel in m ' + str(self.travel) + '\n')
        self.f.write('gradient ' + str(self.gradient) + '\n')
        self.f.write('' + '\n')
        self.f.write('projectile mass in kg ' + str(self.projectile_mass) + '\n')
        self.f.write('switch to calculate if energy lost to air resistance ' + str(self.switch_to_calculate_energy_lost_to_air_resistance) + '\n')
        self.f.write('fraction of work against bore used to heat tube ' + str(self.fraction_of_work_against_bore_to_heat_tube) + '\n')
        self.f.write('gas pressure in front of projectile pa ' + str(self.gas_pressure_in_front_of_projectile) + '\n')
        self.f.write('' + '\n')
        self.f.write('number of barrel resistance points (br,trav) ' + str(self.number_of_barrel_resistance_points) + '\n')
        self.f.write('bore resistance Mpa     travel m' + '\n')
        for resistancePoint in self.barrel_resistance_points:
            self.f.write(' ' + str(resistancePoint.br) + '\t\t\t\t\t\t\t' + str(resistancePoint.trav) + '\n')
        self.f.write('' + '\n')
        self.f.write('mass of recoiling parts kg' + str(self.mass_of_recoiling_parts) + '\n')
        self.f.write('number of recoil points (force,time) should be 2 ' + str(self.number_of_recoiling_parts) + '\n')
        for part in self.recoiling_parts:
            self.f.write(' ' + str(part.force) + '\t\t\t\t\t\t\t' + str(part.time) + '\n')
        self.f.write('' + '\n')
        self.f.write('free convective heat transfer coefficient w/m^2-k ' + str(self.free_convective_heat_transfer_coefficient) + '\n')
        self.f.write('chamber wall thickness m ' + str(self.chamber_wall_thickness) + '\n')
        self.f.write('heat capacity of steel chamber wall j/kg-k ' + str(self.heat_capacity_of_steel_chamber_wall) + '\n')
        self.f.write('initial temperature of chamber wall k ' + str(self.initial_temperature_of_chamber_wall) + '\n')
        self.f.write('heat loss coefficient (should be 1) ' + str(self.heat_loss_coefficient) + '\n')
        self.f.write('density of steel chamber wall kg/m^3 ' + str(self.density_of_steel_chamber_wall) + '\n')
        self.f.write('' + '\n')
        self.f.write('impetus of igniter j/kg ' + str(self.impetus_of_igniter) + '\n')
        self.f.write('covolume of igniter m^3/kg ' + str(self.covolume_of_igniter) + '\n')
        self.f.write('adiabatic flame temperature k ' + str(self.adiabatic_flame_temperature) + '\n')
        self.f.write('mass of igniter kg ' + str(self.mass_of_igniter) + '\n')
        self.f.write('ratio of specific heats for igniter ' + str(self.ratio_of_specific_heats_for_igniter) + '\n')
        self.f.write('' + '\n')
        self.f.write('number of propellants ' + str(self.number_of_propellants) + '\n')
        self.f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            self.f.write('for propellant number ' + str(i) + '\n')
            self.f.write('impetus of the propellant j/kg ' + str(propellant.impetus) + '\n')
            self.f.write('adiabatic flame temperature of propellant k ' + str(propellant.flameTemp) + '\n')
            self.f.write('covolume of the propellant m^3/kg ' + str(propellant.covolume) + '\n')
            self.f.write('mass of propellant kg ' + str(propellant.mass) + '\n')
            self.f.write('density of propellant kg/m^3 ' + str(propellant.density) + '\n')
            self.f.write('ratio of specific heats of propellant ' + str(propellant.ratioOfSpecificHeats) + '\n')
            self.f.write('number of perforations of propellant grain ' + str(propellant.perforations) + '\n')
            self.f.write('length of propellant grain m ' + str(propellant.lengthOfGrain) + '\n')
            self.f.write('diameter of perforation of propellant grain m ' + str(propellant.diameterOfPerforation) + '\n')
            self.f.write('outside diameter of propellant grain m ' + str(propellant.outsideDiameter) + '\n')
        self.f.write('' + '\n')
        i = 0
        for propellant in self.propellants:
            i += 1
            self.f.write(' ' + str(len(propellant.burningRateList)) + ' burning rate points for propellant ' + str(i) + '\n')
            self.f.write('\n')
            self.f.write('exponent   coefficient  pressure' + '\n')
            self.f.write('  -         m/s-MPa^a      MPa' + '\n')
            for burnRate in propellant.burningRateList:
                self.f.write(' ' + str(burnRate.exponent) + '\t\t'+ str(burnRate.coefficient) + '\t\t' + str(burnRate.pressure) + '\n')
        self.f.write('' + '\n')
        self.f.write('time step sec ' + str(self.time_step_sec) + '\n')
        self.f.write('print step sec ' + str(self.print_step_sec) + '\n')
        self.f.write('time to stop (if before projectile exit) sec ' + str(self.time_to_stop) + '\n')

if __name__ == '__main__':
    main('19h')
