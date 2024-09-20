#!/usr/bin/env python

import sys
import pandas as pd

def translateWord(word: str, word_dict: pd.DataFrame, from_col=1, to_col=2):
    trs = word_dict.loc[word_dict[from_col - 1] == word, ][to_col - 1]
    if len(trs) == 0:
        return(None)
    else:
        return(trs.values[0])

def main():
    """
    target_file = 'input/file.tsv'
    dict_file = 'input/dict.tsv'
    from_col = int('1')
    to_col = int('2')
    translate_on = '2,3'
    translate_on = 'all'
    translate_on = '3'
    keep_untranslated = bool(int('1'))
    """
    dict_file = sys.argv[1]
    target_file = sys.argv[2]
    from_col = int(sys.argv[3])
    to_col = int(sys.argv[4])
    translate_on = sys.argv[5]
    keep_untranslated = bool(int(sys.argv[6]))
    
    if translate_on != 'all':
        translate_on = str(translate_on).split(',')
        try:
            translate_on = [int(i) for i in translate_on]
        except:
            translate_on = [int(translate_on)]

    word_dict = pd.read_csv(dict_file, sep='\t', header=None)
    
    if translate_on != 'all':
        with open(target_file, 'r') as target_fh:
            for line in target_fh:
                fields = line.removesuffix('\n').split('\t')
                translated = True
                for i in translate_on:
                    tr_word = translateWord(fields[i-1], word_dict, from_col=from_col, to_col=to_col)
                    if tr_word == None:
                        translated = False
                        break
                    else:
                        fields[i-1] = tr_word
                if translated or keep_untranslated:
                    print('\t'.join([str(field) for field in fields]))
    else:
        with open(target_file, 'r') as target_fh:
            for line in target_fh:
                fields = line.removesuffix('\n').split('\t')
                translated = True
                for i in range(len(fields)):
                    tr_word = translateWord(fields[i-1], word_dict, from_col=from_col, to_col=to_col)
                    if tr_word == None:
                        translated = False
                        break
                    else:
                        fields[i-1] = tr_word
                if translated or keep_untranslated:
                    print('\t'.join([str(field) for field in fields]))  
if __name__ == '__main__':
    main()