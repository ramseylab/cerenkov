"""
Created on Jul 31, 2016

@author: ramseylab
"""

import pandas


def __dict_to_dfm(dictionary, columns, index):
    df = pandas.DataFrame(list(dictionary.items()))
    df.columns = columns
    df = df.set_index(index)
    
    return df


def match(x, x_id, x_pos, y, y_start, y_end, y_value, xy_sorted=False):
    # Initialize a Dictionary for the final result,
    # with `x[x_id]` as keys and empty list as initial values.
    # The final result goes like
            # 'ID_1' : [Value_1]
            # 'ID_2' : [Value_2, Value_3]
    
    result = {key: [] for key in x[x_id]}
    
    if y is None or y.empty:
        return __dict_to_dfm(result, [x_id, y_value], x_id)
    
    if not xy_sorted:
        x.sort_values(by=x_pos, ascending=True, inplace=True)
        y.sort_values(by=y_start, ascending=True, inplace=True)
        
    y_count = y.shape[0]

    # The 'top' index, the first index of `y` that the current `x` matched
    # We must keep this index for the searching of the next x-y match
    i_top = 0
    # The 'bottom' index, the searching index after `iTop` for the current `x`
    # You can also imagine these 2 indices as "top" and "bottom" for a vertical searching window
    i_btm = 0

    # `iterrows` is a generator which yield `(index, row)`
    # I don't need `index` here
    for _, x_row in x.iterrows():
        # Skip all the `y_start` scanned in the previous match
        while i_top < y_count and x_row[x_pos] >= y.iloc[i_top][y_end]:
            i_top += 1

        # No matching until the last `(y_start, y_end)` ranges or
        # the `x_row.x_pos` is in-between two `y` ranges
        if i_top == y_count or (i_top < y_count and x_row[x_pos] < y.iloc[i_top][y_start]):
            # go to next x_row
            continue
        # else:
            # Found a match
        # // end if

        # Memorize the current first matching range index.
        # When matching the next x_row, we'll start from `iTop`,
        # while we'll continue with `i_btm` now
        i_btm = i_top

        # Continue searching the next ranges
        while i_btm < y_count and x_row[x_pos] >= y.iloc[i_btm][y_start]:
            if x_row[x_pos] < y.iloc[i_btm][y_end]:
            # if x_row.y_end <= y.iloc[i_btm].y_end:
                # Found a match
                    # If `extend` is used instead of `append`,
                    # TypeError: 'numpy.int64' object is not iterable
                    result[x_row[x_id]].append(y.iloc[i_btm][y_value])
                # else:
                    # continue
            # else:
                # no match and continue
            # // end if
            i_btm += 1
        # // end while
    # // end for
    
    return __dict_to_dfm(result, [x_id, y_value], x_id)


def group_match(x, x_group_key, x_id, x_pos, y, y_group_key, y_start, y_end, y_value):
    
    def __find_group(group_name, league):
        try:
            return league.get_group(group_name)
        except KeyError:
            return None
    
    x_league = x.groupby(x_group_key)
    y_league = y.groupby(y_group_key)
    
#     for xKeyVal, x_group in xLeague:
#         yGroup = yLeague.get_group(xKeyVal)
#         
#         match(x_group, x_id, x_pos, yGroup, y_start, y_end, y_value, xySorted = False)

    result = x_league.apply(lambda x_group: match(x_group, x_id, x_pos, __find_group(x_group.name, y_league),
                                                  y_start, y_end, y_value, xy_sorted=False))
    result = result.reset_index().set_index(x_id)
    
    return result
