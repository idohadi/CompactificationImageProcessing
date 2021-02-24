function printBegEndMsg(msg, beginningFlag)
if beginningFlag
    fprintf([datestr(datetime), ' ::: Beginning: ', msg, '\n']);
else
    fprintf([datestr(datetime), ' ::: Completed: ', msg, '\n']);

end
end
