#!/usr/bin/env python3
#basic_classes.py

class Circle():
     ''' A simple example of a class with attributes '''
     def __init__(self, color, radius):
          self.color = color
          self.radius = radius
    
     def diameter(self):
          ''' Find the wheel diameter  of the circle '''
          return (self.radius * 2)
          
     def circumference(self):
          ''' Find the circumference of the circle '''
          return (2*3.14*self.radius)
      
     def isRed(self):
          ''' Find the colour of the circle '''
          if self.color == 'red':
             return True
          else:
             return False
class GraduateStudent():
     ''' A simple example of a class with attributes '''
     def __init__(self,first_name,last_name,year,major):
          self.first_name = first_name
          self.last_name = last_name
          self.year = year
          self.major = major
        
     def year_matriculated(self):
          ''' Finds the year a student graduated '''
          if self.year == 1:
             return 2019
          else:
             return (2020 -self.year)

Circle1 = Circle("red", 12)
print(Circle1.color)
print(Circle1.radius)
print(Circle1.diameter())
print(Circle1.circumference())
Circle1.isRed()

graduate_student_1 = GraduateStudent("Mrinal","Subash",6,"Science")
