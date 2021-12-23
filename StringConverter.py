
#{A,T,C,G} = {1,2,3,4]

IndexType = 'P' #M for Matlab, P for Python (1-4 or 0-3 index)

String = 'AGTTGCGT ACGTGCAA'
ConvertedString = ''
KeepSpaces = True
Letters = ['A','T','C','G']

if IndexType == 'M':
    if String[0] in Letters: #We convert letter to index
        for k in range(len(String)):
            if String[k] == 'A':
                ConvertedString += '1'
            elif String[k] == 'T':
                ConvertedString += '2'
            elif String[k] == 'C':
                ConvertedString += '3'
            elif String[k] == 'G':
                ConvertedString += '4'
            elif String[k] == ' ':
                if KeepSpaces:
                    ConvertedString += ' '
                else:
                    ConvertedString += ''
    else:
        for k in range(len(String)):
            if String[k] == '1':
                ConvertedString += 'A'
            elif String[k] == '2':
                ConvertedString += 'T'
            elif String[k] == '3':
                ConvertedString += 'C'
            elif String[k] == '4':
                ConvertedString += 'G'
            elif String[k] == ' ':
                if KeepSpaces:
                    ConvertedString += ' '
                else:
                    ConvertedString += ''

else:
    if String[0] in Letters:  # We convert letter to index
        for k in range(len(String)):
            if String[k] == 'A':
                ConvertedString += '0'
            elif String[k] == 'T':
                ConvertedString += '1'
            elif String[k] == 'C':
                ConvertedString += '2'
            elif String[k] == 'G':
                ConvertedString += '3'
            elif String[k] == ' ':
                if KeepSpaces:
                    ConvertedString += ' '
                else:
                    ConvertedString += ''
    else:
        for k in range(len(String)):
            if String[k] == 'A':
                ConvertedString += '0'
            elif String[k] == 'T':
                ConvertedString += '1'
            elif String[k] == 'C':
                ConvertedString += '2'
            elif String[k] == 'G':
                ConvertedString += '3'
            elif String[k] == ' ':
                if KeepSpaces:
                    ConvertedString += ' '
                else:
                    ConvertedString += ''
print("Original string: " + str(String))
print("Converted String: " + str(ConvertedString))