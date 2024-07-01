import React from "react";
import { TextField, Autocomplete, MenuItem } from "@mui/material";
import CheckIcon from "@mui/icons-material/Check";
import { FixedSizeList } from 'react-window';

const ListboxComponent = React.forwardRef(function ListboxComponent(props, ref) {
    const { children, ...other } = props;
    const itemData = React.Children.toArray(children);
    const ITEM_SIZE = 46;
  
    return (
        <div ref={ref}>
            <FixedSizeList
                height={250}
                width="100%"
                itemSize={ITEM_SIZE}
                itemCount={itemData.length}
                {...other}
                itemData={itemData}
            >
            {({ index, style }) => (
                <div 
                    style={style}
                    // otherwise, the input field will lose focus and you can't select item in list
                    onMouseDown={(event) => event.preventDefault()}
                >
                    {itemData[index]}
                </div>
            )}
            </FixedSizeList>
        </div>
    );
});   

export default function MultiSelect ({ 
    title,
    labels, // as {label: label, value: value} 
    selectedLabels,
    setSelectedLabels
}) {

    const handleChange = (event, value) => {
        setSelectedLabels(value);
    };

    return (
        <Autocomplete
            sx={{ m: 1, width: "100%" }}
            multiple
            options={labels}
            value={selectedLabels}
            onChange={handleChange}
            getOptionLabel={(option) => {
                return option.label.charAt(0).toUpperCase() + option.label.slice(1).replace(/_/g, " ");
            }}
            disableCloseOnSelect
            ListboxComponent={ListboxComponent} // delete if using old render option
            renderInput={(params) => (
                <TextField
                    {...params}
                    variant="outlined"
                    label={title}
                />
            )}
            // renderOption={(props, option, { selected }) => (
            //     <MenuItem
            //         {...props}
            //         key={option.value}
            //         value={option.label}
            //         sx={{ justifyContent: "space-between" }}
            //     >
            //         {option.label.charAt(0).toUpperCase() + option.label.slice(1).replace(/_/g, " ")}
            //         {selected ? <CheckIcon color="info" /> : null}
            //     </MenuItem>
            // )}
            renderOption={(props, option, { selected }) => (
                <MenuItem 
                    {...props} 
                    key={option.value}
                    value={option.label}
                >
                    {option.label.charAt(0).toUpperCase() + option.label.slice(1).replace(/_/g, " ")}
                    {selected ? <CheckIcon color="info" /> : null}
                </MenuItem>
            )}
        />
    );
}