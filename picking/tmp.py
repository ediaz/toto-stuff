# Ask for the number and store it in userNumber
userNumber = input('Give me an integer number: ')

# Make sure the input is an integer number
userNumber = int(userNumber)

# Get the square of the number
# userNumber**2 is the same as saying pow(userNumber, 2)
userNumber = userNumber**2

# Print square of given number
print 'The square of your number is: ' + str(userNumber)

