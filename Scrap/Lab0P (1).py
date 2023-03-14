import random

def main():

    player_ones_num = 1
    player_ones_turn = True
    player_ones_roll = 3

    player_twos_num = 2
    player_twos_turn = False
    player_twos_roll = 3

    if player_ones_turn:
        s = input()
        print("You entered string " + s)

        '''In Java, == compares the memory/reference locations of two things (i.e. variables, Objects)
        If we wish to compare contents (i.e. the values of two Strings), we use .equals()!
        Try using s.equals("/roll") in Java'''
        if s == "/roll":

            '''You can see we are "casting" two ints here as Strings in order to please the Python gods and avoid an error.
            Python is a strongly typed programming language so the + operator can only concatenate Strings.
            In Java, we no longer need to cast ints as Strings in order to concatenate them.'''
            print("Player " + str(player_ones_num) + " Rolling: 0 - " + str(player_twos_roll))

            player_ones_roll = random.randint(0, player_twos_roll)

            print("Player " + str(player_ones_num) + " Rolled A " + str(player_ones_roll))

            player_ones_turn = False
            player_twos_turn = True

    else:
        s = input()
        print("You entered string " + s)

        if s == "/roll":
            print("Player " + str(player_twos_num) + " Rolling: 0 - " + str(player_ones_roll))

            player_twos_roll = random.randint(0, player_ones_roll)

            print("Player " + str(player_twos_num) + " Rolled A " + str(player_twos_roll))

            player_twos_turn = False
            player_ones_turn = True

    if player_ones_roll < player_twos_roll:
        print("Player 1 Wins!")
    elif player_twos_roll < player_ones_roll:
        print("Player 2 Wins!")
    else:
        print("Player 1 And Player 2 Tied!")

main()